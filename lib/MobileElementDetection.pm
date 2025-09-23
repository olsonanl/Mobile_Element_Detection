package MobileElementDetection;

use strict;
use warnings;
use Carp::Always;
use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::Client;
use Bio::KBase::AppService::ReadSet;
use Bio::KBase::AppService::AppSpecs;
use File::Slurp;
use Data::Dumper;
use gjoseqlib;
use JSON::XS;
use IPC::Run;
use POSIX;
use Date::Parse;
use File::Temp;

# Export subroutines in logical order of execution
use Exporter qw(import);
our @EXPORT_OK = qw(
    get_taxonomy_mappings
    preflight
    run
    process_input_and_validate
    validate_fasta_file
    assemble_reads_with_app
    find_app_spec
    await_task_completion
    compute_read_size
    execute_genomad_analysis
    process_viral_contigs
    process_plasmid_contigs
    extract_contig_sequence
    submit_annotation_job
    read_gto_file
    generate_html_report
    organize_results_and_outputs
);
our %EXPORT_TAGS = (
    all => \@EXPORT_OK
);

# Module constants
our $DB_PATH = "/vol/bvbrc/production/application-backend/genomad/genomad_db";
our $ASSEMBLY_APP = "GenomeAssembly2";
our $ANNOTATION_APP = "GenomeAnnotation";

#-----------------------------------------------------------------------------------------
# Taxonomy mappings - keeping this at the top for ease of future updating.
# This is a mapping between the ICTV taxonomy and the corresponding NCBI tax id.
#-----------------------------------------------------------------------------------------

sub get_taxonomy_mappings 
{
    my %lowvan_ids = (
        "Bunyaviricetes"   => 1980410,
        "Paramyxoviridae"  => 11158,
        "Filoviridae"      => 11266,
        "Pneumoviridae"    => 11244,
        "Orthomyxoviridae" => 11308,
    );
    
    my %phano_ids = (
        "Caudoviricetes"   => 2731619,
        "Loebvirae"   => 2732090,
        "Sangervirae"   => 2732091,
        "Vidaverviricetes"   => 2732460,
        "Leviviricetes"   => 2842243,
        "Picobirnaviridae"   => 585893,
        "Partitiviridae"   => 11012,
        "Matsushitaviridae"   => 3152134,
        "Ainoaviricetes"   => 3044425,
        "Vinavirales"   => 2732560,
        "Prepoliviricotina"   => 3412695,
        "Obscuriviridae"   => 3424590,
        "Plasmaviridae"   => 10472,
    );
    
    return (\%lowvan_ids, \%phano_ids);
}

#-----------------------------------------------------------------------------------------
# Resource estimation for job scheduling
#-----------------------------------------------------------------------------------------

sub preflight 
{
    my($app, $app_def, $raw_params, $params) = @_;    
    
    my $cpu = 16;
    my $memory = "32G"; 
    my $runtime;

    if ($params->{input_type} eq 'reads') {
        # Assembly preflight lifted from GenomeAssembly2
        my $readset;
        eval {
            $readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params);
        };
        if ($@) 
        {
            die "Error parsing assembly parameters: $@";
        }

        my($ok, $errs, $comp_size, $uncomp_size) = $readset->validate($app->workspace);

        if (!$ok) 
        {
            die "Reads as defined in parameters failed to validate. Errors:\n\t" . join("\n\t", @$errs) . "\n";
        }
        
        print STDERR "compressed=$comp_size uncompressed=$uncomp_size\n";
		
		#these were lifted from GenomeAssembly2:
        my $est_comp = $comp_size + 0.75 * $uncomp_size;
        $est_comp /= 1e6;  # Convert to megabytes
        
        # Estimated conservative rate is 10sec/MB for compressed data under 1.5G, 4sec/MB for data over that
        my $est_time = int($est_comp < 1500 ? (10 * $est_comp) : (4 * $est_comp));

        # Check for long-read platforms which take longer
        my %plats = map { $_->{platform} => 1 } $readset->libraries;

        if ($plats{nanopore} || $plats{pacbio}) 
        {
            $est_time *= 5;
        }

        # Modern assemblers are more complex - multiply by 10 
        $est_time *= 10;
        $est_time = 3600 if $est_time < 3600;  # Minimum 1 hour

        # Add 4 hours for genomad, this is janky but there is no easy way to estimate 
        # contig number and assembly file size from reads.  This should cover any assembly
        # < 100MB
        $runtime = $est_time + 14400;  #

        # Estimated compressed storage based on input compressed size
        my $est_storage = int(1.3e6 * $est_comp / 0.75);

        $cpu = 32;  # More CPUs for assembly
        $memory = "128G";  # More memory for assembly
        
        # Adjust CPU/memory based on estimated time like GenomeAssembly2
        if ($est_time < 6*60*60) # Less than 6 hours
        {  
            $cpu = 16;
            $memory = "48G";
        }
    }
    else 
    {
        # Contigs: input requirements based on genomad runs
        my $ctg = $params->{input_file};
        $ctg or die "Assembled contigs must be specified\n";
        my $res = $app->workspace->stat($ctg);
        my $size = $res->size;
        $size > 0 or die "Contigs not found\n";
        
        # Runtime estimates based on contigs size
        if ($size < 10_000_000) {
            $runtime = 3600; #1hr wall time for jobs < 10MB
        }
        elsif (($size >= 10_000_000) && ($size < 100_000_000)) {
            $runtime = (3600 * 4); # 4 hr wall time for job < 100MB
        }
        elsif (($size >= 100_000_000) && ($size < 1_000_000_000)) {
            $runtime = (3600 * 12); # 12 hr wall time for job 100MB-1GB
        }
        elsif ($size >= 1_000_000_000) {
            $runtime = (3600 * 72); # 3 day wall time for job > 1GB
        }
    }
    return { cpu => $cpu, memory => $memory, runtime => $runtime };
}

#-----------------------------------------------------------------------------------------
# Main subroutine for mobile element detection workflow
#-----------------------------------------------------------------------------------------

sub run 
{
    my($app, $app_def, $raw_params, $params) = @_;
    
    # Let AppScript handle proper job result folder creation
    my $folder = $app->result_folder();
    
    print "App-Genomad: ", Dumper($app_def, $raw_params, $params);
    my $threads = $ENV{P3_ALLOCATED_CPU} // 16;
    
    # Step 1: Process and validate input
    my ($input_file, $prefix, $contigs_ws_path) = process_input_and_validate($app, $params);
    
    # Step 2: Execute genomad analysis
    execute_genomad_analysis($params, $input_file, $threads);
    
    # Step 3: Read input sequences for processing
    open (IN, "<$input_file") or warn "could not open fasta file of contigs\n";
    my @seqs = &gjoseqlib::read_fasta(\*IN);
    close IN; 
    my $seqH = {};
    for my $i (0..$#seqs) 
    {
        $seqH->{$seqs[$i][0]}->{SEQ} = $seqs[$i][2];
        $seqH->{$seqs[$i][0]}->{ANN} = $seqs[$i][1];
    }
    
    # Step 4: Process viral contigs
    my @viral_data = process_viral_contigs($app, $params, $seqH, $folder, $prefix);
    
    # Step 5: Process plasmid contigs
    my @plasmid_data = process_plasmid_contigs($app, $params, $seqH, $folder, $prefix);
    
    # Step 6: Generate HTML report
    my $html_report = generate_html_report(\@viral_data, \@plasmid_data, $params->{'output_file'}, $folder);
    $app->workspace->save_file_to_file($html_report, {}, "$folder/Analysis_Summary.html", 'html', 1);
    
    # Step 7: Organize results and outputs
    organize_results_and_outputs($app, $params, $folder, $prefix);
}

#-----------------------------------------------------------------------------------------
# Input processing and validation
#-----------------------------------------------------------------------------------------

sub process_input_and_validate 
{
    my($app, $params) = @_;
    
    # Validate input parameters
    if (($params->{input_file}) && ($params->{input_type} eq 'reads')) 
    {
        die "Must declare either contigs or reads\n"
    }
    elsif (($params->{input_type} eq 'reads') && (! $params->{paired_end_libs}) && (! $params->{single_end_libs})&& (! $params->{srr_ids})) 
    {
        die "Must declare reads file(s) or SRA accession(s)\n" 
        #it will likely die upstream of this, but it is here just in case.
    }
    
    my $input_file;
    my $prefix;
    my $contigs_ws_path;
    
    if ($params->{input_type} eq 'reads') 
    {
        # Perform assembly using GenomeAssembly2 app
        $contigs_ws_path = assemble_reads_with_app($app, $params);
        
        # Download the assembled contigs for genomad processing
        $input_file = "assembly_contigs.fasta";
        $app->workspace->download_file($contigs_ws_path, $input_file, 1);
        
        # Extract prefix from the actual filename
        $prefix = $input_file;
        $prefix =~ s/\..+//g;  # Remove extension, so "assembly_contigs.fasta" -> "assembly_contigs"
    }
    else 
    {
        # Use provided contigs
        my $input_path = $params->{'input_file'};
        $input_file = $input_path;
        $input_file =~ s/.+\///g;
        $prefix = $input_file;
        $prefix =~ s/\..+//g; 
        
        $app->workspace->download_file($input_path, $input_file, 1);
    }
    
    # Validate the input file
    validate_fasta_file($input_file);
    
    return ($input_file, $prefix, $contigs_ws_path);
}


#-----------------------------------------------------------------------------------------
# Validate the assembly:  just a quick check to make sure we were given contigs
#-----------------------------------------------------------------------------------------

sub validate_fasta_file 
{
    my($input_file) = @_;
    
    open my $fh, "<", $input_file or die "could not open fasta file of contigs for validation\n";
    my @lines = map { scalar <$fh> } 1..4;
    close $fh;    
    
    die "Input File is not fasta format.\n" unless $lines[0] =~ /^>/; 
    die "Input file is likely FASTQ format. Assembled contigs are required.\n" if grep { /^\+/ } @lines;  
    die "Not DNA FASTA\n" if grep { !/^>/ && /[^ACGTRYKMSWBDHVN\s]/i } @lines;
    
    return 1;
}

#-----------------------------------------------------------------------------------------
# Assembly functions (if reads are provided)
#-----------------------------------------------------------------------------------------

sub assemble_reads_with_app 
{
    my($app, $params) = @_;
    
    print "Starting genome assembly using GenomeAssembly2 app\n";
    
    # Extract assembly-related parameters from the incoming params
    my @assembly_keys = qw(paired_end_libs single_end_libs srr_ids recipe 
                          racon_iter pilon_iter trim target_depth normalize 
                          filtlong genome_size min_contig_len min_contig_cov
                          max_bases debug);
    
    my $assembly_input = { map { exists $params->{$_} ? ($_, $params->{$_}) : () } @assembly_keys };
    
    # Set assembly-specific output parameters
    # Use the job result folder as the base, similar to ComprehensiveGenomeAnalysis.pm
    my $assembly_folder = $app->result_folder . "/assembly";
    $assembly_input->{output_path} = $assembly_folder;
    $assembly_input->{output_file} = "assembly";
    
    $assembly_input->{recipe} ||= 'auto';
    
    print "Assembly parameters: ", Dumper($assembly_input);
    
    my $qtask;
    my $inline = $ENV{P3_CGA_TASKS_INLINE} ? 1 : 0;
    if ($inline) {
        # Run assembly inline
        my $app_spec = find_app_spec($ASSEMBLY_APP);  # Use object-oriented method
        my $tmp = File::Temp->new();
        print $tmp encode_json($assembly_input);
        close($tmp);

        my @cmd = ("App-$ASSEMBLY_APP", "xx", $app_spec, $tmp);
        
        print STDERR "Running inline assembly: @cmd\n";

        my $assembly_start_time = time;
        my $inline_rc = system(@cmd);
        my $assembly_end_time = time;
        if ($inline_rc != 0) 
        {
            die "Inline assembly failed with rc=$inline_rc\n";
        }

        $qtask = 
        {
            id => "assembly_$$",
            app => $ASSEMBLY_APP,
            parameters => $assembly_input,
            user_id => $app->token()->user_id(),
            submit_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $assembly_start_time),
            start_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $assembly_start_time),
            completed_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $assembly_end_time),
            status => 'completed'
        };
    }
    else 
    {
        # Submit assembly as separate job
        my $client = Bio::KBase::AppService::Client->new();
        
        my $task = $client->start_app($ASSEMBLY_APP, $assembly_input, $assembly_folder);
        print "Created assembly task: ", Dumper($task);
        
        my $task_id = $task->{id};
        $qtask = await_task_completion($client, $task_id);
        
        if (!$qtask || $qtask->{status} ne 'completed') {
            die "Assembly failed. Task status: " . ($qtask ? $qtask->{status} : "unknown") . "\n";
        }
    }
    
    print "Assembly completed successfully\n";
    
    # The GenomeAssembly2 service creates the contigs in a .assembly subdirectory
    my $contigs_path = "$assembly_folder/.assembly/assembly_contigs.fasta";
    
    # Verify the contigs were created
    my $stats = eval { $app->workspace->get({ objects => [$contigs_path], metadata_only => 1 }) };
    if (!$stats || @$stats == 0) 
    {
        die "Could not find generated contigs in $contigs_path\n";
    }
    
    print "Assembly contigs available at: $contigs_path\n";
    return $contigs_path;
}


#-----------------------------------------------------------------------------------------
# Find the App Spec (GenomeAssembly and GenomeAnnotion need to be checked out from here:
#  git clone https://github.com/bv-brc/p3_assembly
#  git clone https://github.com/bv-brc/p3_genome_annotation
#-----------------------------------------------------------------------------------------

sub find_app_spec 
{
    my($app_name) = @_;
    my $specs = Bio::KBase::AppService::AppSpecs->new;
    my($spec, $spec_file) = $specs->find($app_name);
    return $spec_file;
}

#-----------------------------------------------------------------------------------------
# Await task completion
#-----------------------------------------------------------------------------------------

sub await_task_completion 
{
    my($client, $task_id, $query_frequency, $timeout) = @_;

    $query_frequency //= 30;  # Check every 30 seconds

    my %final_states = map { $_ => 1 } qw(failed suspend completed user_skipped skipped passed deleted);

    my $end_time;
    if ($timeout) 
    {
        $end_time = time + $timeout;
    }

    my $qtask;
    print "Waiting for task $task_id to complete...\n";
    
    while (!$end_time || (time < $end_time)) 
    {
        my $qtasks = eval { $client->query_tasks([$task_id]); };
        if ($@) 
        {
            warn "Error checking task status: $@\n";
        }
        else 
        {
            $qtask = $qtasks->{$task_id};
            my $status = $qtask->{status};
            print "Task $task_id status: $status\n";
            
            last if $final_states{$status};
        }
        sleep($query_frequency);
        undef $qtask;
    }
        
    if (!$qtask) {
        die "Task $task_id did not complete within timeout\n";
    }
    
    return $qtask;
}

#-----------------------------------------------------------------------------------------
# Compute read set size
# Lifted from GenomeAssembly2
#-----------------------------------------------------------------------------------------

sub compute_read_size 
{
    my($app, $params) = @_;
    
    # Use ReadSet module like GenomeAssembly2 does
    my $readset;
    eval 
    {
        $readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params);
    };
    if ($@) 
    {
        warn "Error creating ReadSet: $@";
        return 0;
    }
    
    my($ok, $errs, $comp_size, $uncomp_size) = $readset->validate($app->workspace);
    
    if (!$ok) 
    {
        warn "ReadSet validation failed: " . join(", ", @$errs);
        return 0;
    }
    
    # Use the estimated compressed size + 75% of uncompressed size
    # This was lifted from from GenomeAssembly2
    my $estimated_size = $comp_size + 0.75 * $uncomp_size;
    
    print STDERR "ReadSet validation: comp_size=$comp_size, uncomp_size=$uncomp_size, estimated=$estimated_size\n";
    
    return $estimated_size;
}

#-----------------------------------------------------------------------------------------
# GeNomad execution
#-----------------------------------------------------------------------------------------

sub execute_genomad_analysis 
{
    my($params, $input_file, $threads) = @_;
    
    my @command;
    my $of = $params->{'output_file'};
    
    push @command, ("genomad", "end-to-end", "--threads", "$threads"); 
    if ($params->{'cleanup'}){push @command, "--cleanup";}
    if ($params->{'restart'}){push @command, "--restart";}
    if ($params->{'verbose'}){push @command, "--verbose";}
    if ($params->{'lenient-taxonomy'}){push @command, "--lenient-taxonomy";}
    if ($params->{'full-ictv-lineage'}){push @command, "--full-ictv-lineage";}
    if ($params->{'force-auto'}){push @command, "--force-auto";}
    if ($params->{'filtering-preset'}){push @command, "--$params->{'filtering-preset'}";}
    if ($params->{'composition'}) 
    {    
        push @command, "--composition";
        push @command, $params->{'composition'};
    }
    push @command, $input_file;
    push @command, $of;
    push @command, $DB_PATH;
    
    my $genomad_rc = system(@command);   
    die "Failure running genomad" if $genomad_rc != 0;
    
    return 1;
}

#-----------------------------------------------------------------------------------------
# Results processing - viral contigs
#-----------------------------------------------------------------------------------------

sub process_viral_contigs 
{
    my($app, $params, $seqH, $folder, $prefix) = @_;
    
    my $of = $params->{'output_file'};
    my @viral_data = ();
    my ($lowvan_ids, $phano_ids) = get_taxonomy_mappings();
    
    return @viral_data unless (-s "$of/${prefix}_summary/${prefix}_virus_summary.tsv");
    
    open (IN, "<$of/${prefix}_summary/${prefix}_virus_summary.tsv") or warn "could not open viral TSV file\n"; 
    
    while (<IN>) 
    {
        chomp;
        next if /^seq_name/; # Skip header line
        my ($id, $len, $topo, $coords, $genes, $code, $score, $fdr, $hallmarks, $enrichment, $tax) = split /\t/; 
        my @labels = split (";", $tax);
        
        my $annotated = undef;
        my $recipe = undef;
        my $tax_id = undef;
        
        # Process each viral contig that meets the LowVan or Phannotate criteria        
        foreach my $label (@labels) 
        {
            if ((exists $lowvan_ids->{$label}) || (exists $phano_ids->{$label})) 
            {
                # Determine annotation type and recipe
                if (exists $lowvan_ids->{$label}) 
                {
                    $annotated = "LowVan"; 
                    $recipe = "viral-lowvan";
                    $tax_id = $lowvan_ids->{$label};
                }
                elsif (exists $phano_ids->{$label}) 
                {
                    $annotated = "Phannotate"; 
                    $recipe = "phage";
                    $tax_id = $phano_ids->{$label};
                }

                # Create file name
                my $fn = $label . "_" . $id;
                $fn =~ s/[^A-Za-z0-9]/_/g;    
                
                # Extract sequence (handles both regular contigs and proviruses)
                my ($contig_to_use, $dna_string, $contig_len) = extract_contig_sequence($seqH, $id);
                
                # Set up annotation parameters
                my $annotation_params = 
                {
                    output_dir => "$of/Viral_Annotation",
                    workspace_base_dir => "$folder/Viral_Annotation",
                    code => 11,
                    domain => "Viruses",
                    scientific_name => $label,
                    recipe => $recipe,
                    taxonomy_id => $tax_id,
                };
                
                # Submit annotation job
                my $annotation_success = submit_annotation_job($app, $folder, $fn, $id, $dna_string, $annotation_params);
                
                if ($annotation_success) 
                {
                    # Download and process GTO file
                    unless (-e "$of/Viral_GTO") { mkdir "$of/Viral_GTO"; }
                    my $gto = "$folder/Viral_Annotation/$fn/\.$fn/$fn.genome";
                    my $local = "$of/Viral_GTO/$fn.genome";                
                    $app->workspace->download_file($gto, $local, 1);
                    
                    # Read GTO file and collect data for HTML report
                    my $gto_data = read_gto_file($local);
                    $tax =~ s/;+$//;
                    $tax =~ s/^Viruses;//g;
                    
                    my $viral_row = 
                    {
                        contig => $contig_to_use,
                        file_name => $fn,
                        anno => $annotated,
                        genome_id => $gto_data ? $gto_data->{id} : 'N/A',
                        length => $contig_len,
                        gc_content => $gto_data ? $gto_data->{gc_content} : 'N/A',
                        genome_quality => ($annotated eq "LowVan" && $gto_data) ? $gto_data->{genome_quality} : 'N/A',
                        topology => $topo,
                        score => $score,
                        bvbrc_genes => $gto_data ? $gto_data->{feature_count} : 'N/A',
                        taxonomy => $tax
                    };
                    
                    push @viral_data, $viral_row;
                }
                last; # Only process the first matching taxonomy label per contig
            }                
        }
        
        # If no LowVan or Phannotate annotation was done, still add to viral data
        if (!$annotated) 
        {
            my $viral_row = 
            {
                contig => $id,
                file_name => 'N/A',
                anno => 'N/A',
                genome_id => 'N/A',
                length => $len,
                gc_content => 'N/A',
                genome_quality => 'N/A',
                topology => $topo,
                score => $score,
                bvbrc_genes => 'N/A',
                taxonomy => $tax
            };               
            push @viral_data, $viral_row;
        }
    }
    close(IN);
    
    unless (-e "$of/Viral_Annotation") 
    {
        print STDERR "No Phage or LowVan supported viral taxa identified\n";
    }
    
    return @viral_data;
}


#-----------------------------------------------------------------------------------------
# Results processing - plasmid contigs
#-----------------------------------------------------------------------------------------

sub process_plasmid_contigs 
{
    my($app, $params, $seqH, $folder, $prefix) = @_;
    
    my $of = $params->{'output_file'};
    my @plasmid_data = ();
    
    return @plasmid_data unless (-s "$of/${prefix}_summary/${prefix}_plasmid_summary.tsv");
    
    open (IN, "<$of/${prefix}_summary/${prefix}_plasmid_summary.tsv") or warn "could not open plasmid TSV file\n"; 
    
    while (<IN>) 
    {
        chomp;
        next if /^seq_name/; # Skip header line
        my ($id, $len, $topo, $genes, $code, $score, $fdr, $hallmarks, $enrichment, $cnj_genes, $amr_genes) = split /\t/; 
        
        next unless $len >= 300; # Hardcoded min len should match min contig len in the auto assembly recipe
        
        my $annotated = "RAST";
        my $fn = $id; # set file name equal to the contig id
        $fn =~ s/[^A-Za-z0-9]/_/g;    
        
        # Set up annotation parameters
        my $annotation_params = {
            output_dir => "$of/Plasmid_Annotation",
            workspace_base_dir => "$folder/Plasmid_Annotation",
            code => 11,
            domain => "Bacteria",
            scientific_name => "Bacterial plasmid",
            recipe => "default",
            taxonomy_id => "2",
        };
        
        # Submit annotation job
        my $annotation_success = submit_annotation_job($app, $folder, $fn, $id, $seqH->{$id}->{SEQ}, $annotation_params);
        
        if ($annotation_success) 
        {
            # Download and process GTO file
            unless (-e "$of/Plasmid_GTO") { mkdir "$of/Plasmid_GTO"; }
            my $gto = "$folder/Plasmid_Annotation/$fn/\.$fn/$fn.genome";
            my $local = "$of/Plasmid_GTO/$fn.genome";                
            $app->workspace->download_file($gto, $local, 1);
            
            # Read GTO file and collect data for HTML report
            my $gto_data = read_gto_file($local);
            
            my $plasmid_row = 
            {
                contig => $id,
                file_name => $fn,
                anno => $annotated,
                genome_id => $gto_data ? $gto_data->{id} : 'N/A',
                length => $len,
                gc_content => $gto_data ? $gto_data->{gc_content} : 'N/A',
                topology => $topo,
                score => $score,
                bvbrc_genes => $gto_data ? $gto_data->{feature_count} : 'N/A',
                taxonomy => 'N/A'  # Plasmids don't have ICTV taxonomy
            };
            push @plasmid_data, $plasmid_row;
        }
    }
    close(IN);
    
    return @plasmid_data;
}

#-----------------------------------------------------------------------------------------
# Extract contig sequence
#-----------------------------------------------------------------------------------------

sub extract_contig_sequence {
    my($seqH, $id) = @_;
    
    if ($id =~ /provirus/) {
        # Example provirus id: accn|BA000007|provirus_1245607_1309390
        my ($contig_id, $prov_start, $prov_end) = $id =~ /(.*)\|provirus_(\d+)_(\d+)$/;
        my $offset = $prov_start - 1;
        my $prov_len = $prov_end - $prov_start + 1;
        my $dna_string = substr($seqH->{$contig_id}->{SEQ}, $offset, $prov_len);
        
        return ($contig_id, $dna_string, $prov_len);
    }
    else 
    {
        return ($id, $seqH->{$id}->{SEQ}, length($seqH->{$id}->{SEQ}));
    }
}


#-----------------------------------------------------------------------------------------
# Submit annotation
#-----------------------------------------------------------------------------------------

sub submit_annotation_job 
{
    my($app, $folder, $fn, $id, $sequence, $annotation_params) = @_;
    
    # Create the contig file
    my $output_dir = $annotation_params->{output_dir};
    unless (-e $output_dir) { mkdir $output_dir; }
    mkdir "$output_dir/$fn";
    
    open (OUT, ">$output_dir/$fn/$fn.contigs.fasta");
    &gjoseqlib::print_alignment_as_fasta(\*OUT, ([$id, "", $sequence]));
    close OUT;
    
    # Create workspace directories
    eval { 
        $app->workspace->create({
            overwrite => 1, 
            objects => [[$annotation_params->{workspace_base_dir}, 'folder']]
        }); 
    };
    eval { 
        $app->workspace->create({
            overwrite => 1, 
            objects => [["$annotation_params->{workspace_base_dir}/$fn", 'folder']]
        }); 
    };
        
    $app->workspace->save_file_to_file("$output_dir/$fn/$fn.contigs.fasta", {}, "$annotation_params->{workspace_base_dir}/$fn/$fn.contigs.fasta", 'contigs', 1);

    # Set up the annotation parameters
    my $descr = 
    {
        contigs => "$annotation_params->{workspace_base_dir}/$fn/$fn.contigs.fasta",
        code => $annotation_params->{code},
        domain => $annotation_params->{domain},
        scientific_name => $annotation_params->{scientific_name},
        output_path => "$annotation_params->{workspace_base_dir}/$fn",
        output_file => $fn,
        queue_nowait => 1,
        analyze_quality => 0,
        skip_indexing => 0,
        recipe => $annotation_params->{recipe},
        taxonomy_id => $annotation_params->{taxonomy_id},
        enable_debug => 0,
    };
    
    # Add recipe-specific parameters
    if ($annotation_params->{recipe} eq "viral-lowvan") 
    {
        $descr->{lowvan_min_contig_length} = 300;
        $descr->{lowvan_max_contig_length} = 30000;
    }
        
    my $json = JSON::XS->new->pretty->canonical();            
    my $tmp = File::Temp->new;
    print $tmp $json->encode($descr);
    close($tmp);
        
    # Submit annotation job
    my $anno_spec = find_app_spec($ANNOTATION_APP);
    my @cmd = ("App-$ANNOTATION_APP", "xx", $anno_spec, "$tmp");
    print STDERR "Run annotation for $fn: @cmd\n";
    
    my $annotation_start_time = time;
    my $anno_rc = system(@cmd);
    my $annotation_end_time = time;
    
    if ($anno_rc != 0) 
    {
        warn "Annotation failed for $fn with rc=$anno_rc\n";
        return undef;
    }
    
    print STDERR "Annotation completed for $fn\n";
    return 1;
}

#-----------------------------------------------------------------------------------------
# Read Annotation GTO file
#-----------------------------------------------------------------------------------------

sub read_gto_file 
{
    my($gto_path) = @_;
    
    unless (-e $gto_path) 
    {
        warn "GTO file not found: $gto_path\n";
        return undef;
    }
        
    my $gto_content = eval { read_file($gto_path) };
    if ($@) 
    {
        warn "Error reading GTO file $gto_path: $@\n";
        return undef;
    }
    
    my $gto = eval { decode_json($gto_content) };
    if ($@) 
    {
        warn "Error parsing JSON in GTO file $gto_path: $@\n";
        return undef;
    }
    
    my $result = 
    {
        id => $gto->{id} || "N/A",
        scientific_name => $gto->{scientific_name} || "N/A",
        gc_content => $gto->{quality}->{gc_content} || "N/A",
        genome_quality => $gto->{quality}->{genome_quality} || "N/A",
        feature_count => 0
    };
    
    # Count features 
    if ($gto->{features} && ref($gto->{features}) eq 'ARRAY') 
    {
        $result->{feature_count} = scalar(@{$gto->{features}});
    }
    return $result;
}

#-----------------------------------------------------------------------------------------
# HTML report generation
#-----------------------------------------------------------------------------------------

sub generate_html_report 
{
    my($viral_data, $plasmid_data, $output_dir, $workspace_folder) = @_;
    
    my $html_file = "$output_dir/mobile_element_detection_report.html";
    
    open(my $html_fh, '>', $html_file) or warn "Cannot open $html_file for writing: $!";
    
    print $html_fh <<'EOF';
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Mobile Element Detection Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            color: #2c3e50;
            text-align: center;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 30px;
            margin-bottom: 15px;
            border-left: 4px solid #3498db;
            padding-left: 10px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 30px;
            font-size: 12px;
            table-layout: fixed;
        }
        th {
            background-color: #3498db;
            color: white;
            padding: 12px 8px;
            text-align: left;
            font-weight: bold;
            cursor: pointer;
            position: relative;
            user-select: none;
            word-wrap: break-word;
            overflow-wrap: break-word;
            border-right: 1px solid #2980b9;
        }
        th:hover {
            background-color: #2980b9;
        }
        th::after {
            content: ' ↕';
            font-size: 10px;
            opacity: 0.5;
            margin-left: 5px;
        }
        th.sort-asc::after {
            content: ' ↑';
            opacity: 1;
        }
        th.sort-desc::after {
            content: ' ↓';
            opacity: 1;
        }
        .resize-handle {
            position: absolute;
            top: 0;
            right: 0;
            width: 5px;
            height: 100%;
            cursor: col-resize;
            background: transparent;
            z-index: 1;
        }
        .resize-handle:hover {
            background-color: rgba(255,255,255,0.3);
        }
        td {
            padding: 10px 8px;
            border-bottom: 1px solid #ddd;
            word-wrap: break-word;
            overflow-wrap: break-word;
            max-width: 0;
            border-right: 1px solid #eee;
        }
        
        /* Set initial column widths for viral table */
        #viralTable th:nth-child(1), #viralTable td:nth-child(1) { width: 9%; }  /* Contig */
        #viralTable th:nth-child(2), #viralTable td:nth-child(2) { width: 13%; } /* File Name */
        #viralTable th:nth-child(3), #viralTable td:nth-child(3) { width: 7%; }  /* Annotation */
        #viralTable th:nth-child(4), #viralTable td:nth-child(4) { width: 9%; }  /* Genome ID */
        #viralTable th:nth-child(5), #viralTable td:nth-child(5) { width: 7%; }  /* Length */
        #viralTable th:nth-child(6), #viralTable td:nth-child(6) { width: 7%; }  /* GC Content */
        #viralTable th:nth-child(7), #viralTable td:nth-child(7) { width: 9%; }  /* Genome Quality */
        #viralTable th:nth-child(8), #viralTable td:nth-child(8) { width: 7%; }  /* Topology */
        #viralTable th:nth-child(9), #viralTable td:nth-child(9) { width: 7%; }  /* Score */
        #viralTable th:nth-child(10), #viralTable td:nth-child(10) { width: 7%; } /* BVBRC */
        #viralTable th:nth-child(11), #viralTable td:nth-child(11) { width: 35%; } /* Taxonomy */
        
        /* Set initial column widths for plasmid table */
        #plasmidTable th:nth-child(1), #plasmidTable td:nth-child(1) { width: 12%; } /* Contig */
        #plasmidTable th:nth-child(2), #plasmidTable td:nth-child(2) { width: 16%; } /* File Name */
        #plasmidTable th:nth-child(3), #plasmidTable td:nth-child(3) { width: 10%; } /* Annotation */
        #plasmidTable th:nth-child(4), #plasmidTable td:nth-child(4) { width: 14%; } /* Genome ID */
        #plasmidTable th:nth-child(5), #plasmidTable td:nth-child(5) { width: 12%; } /* Length */
        #plasmidTable th:nth-child(6), #plasmidTable td:nth-child(6) { width: 12%; } /* GC Content */
        #plasmidTable th:nth-child(7), #plasmidTable td:nth-child(7) { width: 10%; } /* Topology */
        #plasmidTable th:nth-child(8), #plasmidTable td:nth-child(8) { width: 10%; } /* Score */
        #plasmidTable th:nth-child(9), #plasmidTable td:nth-child(9) { width: 14%; } /* BVBRC */
        
        tr:nth-child(even) {
            background-color: #f9f9f9;
        }
        tr:hover {
            background-color: #e8f4f8;
        }
        .genome-link, .feature-link {
            color: #3498db;
            text-decoration: none;
            font-weight: bold;
        }
        .genome-link:hover, .feature-link:hover {
            text-decoration: underline;
        }
        .no-data {
            text-align: center;
            font-style: italic;
            color: #7f8c8d;
            padding: 20px;
        }
        .summary {
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 20px;
        }
        .resize-info {
            font-size: 11px;
            color: #7f8c8d;
            margin-bottom: 10px;
            font-style: italic;
        }
    </style>
    <script>
        let isResizing = false;
        let currentTable = null;
        let currentColumn = null;
        let startX = 0;
        let startWidth = 0;

        function initializeResizableColumns() {
            const tables = document.querySelectorAll('table');
            tables.forEach(table => {
                const headers = table.querySelectorAll('th');
                headers.forEach((header, index) => {
                    // Don't add resize handle to last column
                    if (index < headers.length - 1) {
                        const resizeHandle = document.createElement('div');
                        resizeHandle.className = 'resize-handle';
                        resizeHandle.addEventListener('mousedown', initResize);
                        header.appendChild(resizeHandle);
                    }
                });
            });
        }

        function initResize(e) {
            isResizing = true;
            const resizeHandle = e.target;
            const header = resizeHandle.parentElement;
            currentTable = header.closest('table');
            currentColumn = Array.from(header.parentElement.children).indexOf(header);
            startX = e.pageX;
            startWidth = parseInt(document.defaultView.getComputedStyle(header).width, 10);
            
            document.addEventListener('mousemove', doResize);
            document.addEventListener('mouseup', stopResize);
            e.preventDefault();
        }

        function doResize(e) {
            if (!isResizing) return;
            
            const width = (startWidth + e.pageX - startX);
            if (width > 50) { // Minimum width of 50px
                const tableWidth = currentTable.offsetWidth;
                const widthPercent = (width / tableWidth) * 100;
                
                // Update the width of all cells in this column
                const rows = currentTable.querySelectorAll('tr');
                rows.forEach(row => {
                    const cell = row.children[currentColumn];
                    if (cell) {
                        cell.style.width = widthPercent + '%';
                    }
                });
            }
        }

        function stopResize() {
            isResizing = false;
            currentTable = null;
            currentColumn = null;
            document.removeEventListener('mousemove', doResize);
            document.removeEventListener('mouseup', stopResize);
        }

        function sortTable(tableId, columnIndex) {
            const table = document.getElementById(tableId);
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr'));
            const header = table.querySelector('thead').querySelectorAll('th')[columnIndex];
            
            // Remove sort classes from all headers
            table.querySelectorAll('th').forEach(th => {
                th.classList.remove('sort-asc', 'sort-desc');
            });
            
            // Determine sort direction
            const currentSort = header.getAttribute('data-sort') || 'none';
            const isAsc = currentSort !== 'asc';
            
            // Sort rows
            rows.sort((a, b) => {
                const aVal = a.cells[columnIndex].textContent.trim();
                const bVal = b.cells[columnIndex].textContent.trim();
                
                // Try to parse as numbers first
                const aNum = parseFloat(aVal.replace(/[^0-9.-]/g, ''));
                const bNum = parseFloat(bVal.replace(/[^0-9.-]/g, ''));
                
                let comparison = 0;
                if (!isNaN(aNum) && !isNaN(bNum)) {
                    comparison = aNum - bNum;
                } else {
                    comparison = aVal.localeCompare(bVal);
                }
                
                return isAsc ? comparison : -comparison;
            });
            
            // Update header sort indicator
            header.classList.add(isAsc ? 'sort-asc' : 'sort-desc');
            header.setAttribute('data-sort', isAsc ? 'asc' : 'desc');
            
            // Reorder rows in DOM
            rows.forEach(row => tbody.appendChild(row));
        }

        // Initialize column resizing when page loads
        document.addEventListener('DOMContentLoaded', initializeResizableColumns);
    </script>
</head>
<body>
    <div class="container">
        <h1>Mobile Element Detection Report</h1>
        
        <div class="summary">
            <strong>Summary:</strong> 
EOF
    
    my $viral_count = scalar(@$viral_data);
    my $plasmid_count = scalar(@$plasmid_data);
    print $html_fh "$viral_count viral contigs and $plasmid_count plasmid contigs were found.\n";
    print $html_fh "        </div>\n\n";
    
    # Viral table
    print $html_fh "        <h2>Viral Contigs ($viral_count)</h2>\n";
    if (@$viral_data) 
    {
        print $html_fh <<'EOF';
        <table id="viralTable">
            <thead>
                <tr>
                    <th onclick="sortTable('viralTable', 0)">Contig</th>
                    <th onclick="sortTable('viralTable', 1)">File Name</th>
                    <th onclick="sortTable('viralTable', 2)">Annotation Type</th>
                    <th onclick="sortTable('viralTable', 3)">Genome ID</th>
                    <th onclick="sortTable('viralTable', 4)">Contig Length</th>
                    <th onclick="sortTable('viralTable', 5)">GC Content (%)</th>
                    <th onclick="sortTable('viralTable', 6)">Genome Quality</th>
                    <th onclick="sortTable('viralTable', 7)">Topology</th>
                    <th onclick="sortTable('viralTable', 8)">GeNomad Score</th>
                    <th onclick="sortTable('viralTable', 9)">BVBRC Gene Count</th>
                    <th onclick="sortTable('viralTable', 10)">ICTV Taxonomy</th>
                </tr>
            </thead>
            <tbody>
EOF
        
        foreach my $row (@$viral_data) 
        {
            # Use regular text for annotation type instead of badges
            my $annotation_text = $row->{anno} || 'N/A';
            
            my $genome_link = $row->{genome_id} ne 'N/A' ? 
                "<a href=\"https://www.bv-brc.org/view/Genome/$row->{genome_id}\" class=\"genome-link\" target=\"_blank\">$row->{genome_id}</a>" : 
                'N/A';
            
            my $gc_content = $row->{gc_content} ne 'N/A' ? 
                sprintf("%.2f", $row->{gc_content}) : 'N/A';
            
            # Create hyperlinked BVBRC gene count
            my $bvbrc_genes_link = ($row->{genome_id} ne 'N/A' && $row->{bvbrc_genes} ne 'N/A') ?
                "<a href=\"https://www.bv-brc.org/view/Genome/$row->{genome_id}#view_tab=features\" class=\"feature-link\" target=\"_blank\">$row->{bvbrc_genes}</a>" :
                $row->{bvbrc_genes};
            
            # Create hyperlinked file name for viral annotation directory
            my $file_name_display = $row->{file_name};
            if ($row->{file_name} ne 'N/A' && $workspace_folder) 
            {
                my $annotation_url = "https://www.bv-brc.org/workspace$workspace_folder/Viral_Annotation/$row->{file_name}/$row->{file_name}";
                $file_name_display = "<a href=\"$annotation_url\" class=\"genome-link\" target=\"_blank\">$row->{file_name}</a>";
            }
            
            print $html_fh "                <tr>\n";
            print $html_fh "                    <td>$row->{contig}</td>\n";
            print $html_fh "                    <td>$file_name_display</td>\n";
            print $html_fh "                    <td>$annotation_text</td>\n";
            print $html_fh "                    <td>$genome_link</td>\n";
            print $html_fh "                    <td>$row->{length}</td>\n";
            print $html_fh "                    <td>$gc_content</td>\n";
            print $html_fh "                    <td>$row->{genome_quality}</td>\n";
            print $html_fh "                    <td>$row->{topology}</td>\n";
            print $html_fh "                    <td>$row->{score}</td>\n";
            print $html_fh "                    <td>$bvbrc_genes_link</td>\n";
            print $html_fh "                    <td>$row->{taxonomy}</td>\n";
            print $html_fh "                </tr>\n";
        }
        
        print $html_fh "            </tbody>\n        </table>\n";
    } 
    else 
    {
        print $html_fh "        <div class=\"no-data\">No viral contigs found.</div>\n";
    }
    
    # Plasmid table (no annotation type or ICTV taxonomy columns)
    print $html_fh "\n        <h2>Plasmid Contigs ($plasmid_count)</h2>\n";
    if (@$plasmid_data) 
    {
        print $html_fh <<'EOF';
        <table id="plasmidTable">
            <thead>
                <tr>
                    <th onclick="sortTable('plasmidTable', 0)">Contig</th>
                    <th onclick="sortTable('plasmidTable', 1)">File Name</th>
                    <th onclick="sortTable('plasmidTable', 2)">Annotation Type</th>
                    <th onclick="sortTable('plasmidTable', 3)">Genome ID</th>
                    <th onclick="sortTable('plasmidTable', 4)">Contig Length</th>
                    <th onclick="sortTable('plasmidTable', 5)">GC Content (%)</th>
                    <th onclick="sortTable('plasmidTable', 6)">Topology</th>
                    <th onclick="sortTable('plasmidTable', 7)">GeNomad Score</th>
                    <th onclick="sortTable('plasmidTable', 8)">BVBRC Gene Count</th>
                </tr>
            </thead>
            <tbody>
EOF
        
        foreach my $row (@$plasmid_data) 
        {
            # Use regular text for annotation type
            my $annotation_text = $row->{anno} || 'N/A';
            
            my $genome_link = $row->{genome_id} ne 'N/A' ? 
                "<a href=\"https://www.bv-brc.org/view/Genome/$row->{genome_id}\" class=\"genome-link\" target=\"_blank\">$row->{genome_id}</a>" : 
                'N/A';
            
            my $gc_content = $row->{gc_content} ne 'N/A' ? 
                sprintf("%.2f", $row->{gc_content}) : 'N/A';
            
            # Create hyperlinked BVBRC gene count
            my $bvbrc_genes_link = ($row->{genome_id} ne 'N/A' && $row->{bvbrc_genes} ne 'N/A') ?
                "<a href=\"https://www.bv-brc.org/view/Genome/$row->{genome_id}#view_tab=features\" class=\"feature-link\" target=\"_blank\">$row->{bvbrc_genes}</a>" :
                $row->{bvbrc_genes};
            
            # Create hyperlinked file name for plasmid annotation directory
            my $file_name_display = $row->{file_name};
            if ($row->{file_name} ne 'N/A' && $workspace_folder) 
            {
                my $annotation_url = "https://www.bv-brc.org/workspace$workspace_folder/Plasmid_Annotation/$row->{file_name}/$row->{file_name}";
                $file_name_display = "<a href=\"$annotation_url\" class=\"genome-link\" target=\"_blank\">$row->{file_name}</a>";
            }
            
            print $html_fh "                <tr>\n";
            print $html_fh "                    <td>$row->{contig}</td>\n";
            print $html_fh "                    <td>$file_name_display</td>\n";
            print $html_fh "                    <td>$annotation_text</td>\n";
            print $html_fh "                    <td>$genome_link</td>\n";
            print $html_fh "                    <td>$row->{length}</td>\n";
            print $html_fh "                    <td>$gc_content</td>\n";
            print $html_fh "                    <td>$row->{topology}</td>\n";
            print $html_fh "                    <td>$row->{score}</td>\n";
            print $html_fh "                    <td>$bvbrc_genes_link</td>\n";
            print $html_fh "                </tr>\n";
        }
        
        print $html_fh "            </tbody>\n        </table>\n";
    } 
    else 
    {
        print $html_fh "        <div class=\"no-data\">No plasmid contigs found.</div>\n";
    }
    
    print $html_fh <<'EOF';
    </div>
</body>
</html>
EOF
    
    close($html_fh);
    print STDERR "HTML report generated: $html_file\n";
    return $html_file;
}

#-----------------------------------------------------------------------------------------
# Output organization and cleanup
#-----------------------------------------------------------------------------------------

sub organize_results_and_outputs {
    my($app, $params, $folder, $prefix) = @_;
    
    my $of = $params->{'output_file'};
    
    # Concatenate log files
    my $cat = "cat $of/${prefix}_annotate.log $of/${prefix}_find_proviruses.log $of/${prefix}_marker_classification.log $of/${prefix}_nn_classification.log $of/${prefix}_aggregated_classification.log $of/${prefix}_summary.log > $of/geNomad_run.stderr";

    my $cat_rc = system($cat);   
    die "Failure concatenating log file" if $cat_rc != 0;

    # Create prokka directory for reorganization
    eval { 
        $app->workspace->create({
        overwrite => 1, 
        objects => [["$folder/prokka", 'folder']]
        }); 
    };
    
    # Reorganize output files
    my @to_move = opendir (DIR, "$of/${prefix}_summary/");
    my @reorg = grep{$_ !~ /^\./ && $_ !~ /log/ }readdir(DIR); 
    closedir DIR;
    
    foreach (@reorg) {
        if ($_ =~ /summary\.tsv$/) {
            $app->workspace->save_file_to_file("$of/${prefix}_summary/$_", {}, "$folder/$_", 'tsv', 1);
        }
        elsif ($_ =~ /genes\.tsv$/) {
            $app->workspace->save_file_to_file("$of/${prefix}_summary/$_", {}, "$folder/prokka/$_", 'tsv', 1);
        }
        elsif ($_ =~ /\.faa$/) {
            $app->workspace->save_file_to_file("$of/${prefix}_summary/$_", {}, "$folder/prokka/$_", 'feature_protein_fasta', 1);
        }
        elsif ($_ =~ /\.fna$/) {
            $app->workspace->save_file_to_file("$of/${prefix}_summary/$_", {}, "$folder/$_", 'contigs', 1);
        }
        elsif ($_ =~ /\.json$/) {
            $app->workspace->save_file_to_file("$of/${prefix}_summary/$_", {}, "$folder/$_", 'json', 1);
        }
    }
    
    $app->workspace->save_file_to_file("$of/geNomad_run.stderr", {}, "$folder/geNomad_run.stderr", 'txt', 1);
}
1;  # Required for Perl modules