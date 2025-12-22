import os
import re
import sys
from collections import defaultdict


locus_name_map = {
    "gb|J01866.1|HUMRRB Human 5.8S ribosomal RNA": "rRNA_5.8S",
    "ref|NR_003287.4| Homo sapiens RNA, 28S ribosomal N5 (RNA28SN5), ribosomal RNA": "rRNA_28S",
    "ref|NR_023363.1| Homo sapiens RNA, 5S ribosomal 1 (RNA5S1), ribosomal RNA": "rRNA_5S",
    "ref|NR_046235.3| Homo sapiens RNA, 45S pre-ribosomal N5 (RNA45SN5), ribosomal RNA": "rRNA_45S",
    "ref|NR_137294.1| Homo sapiens mitochondrially encoded 12S ribosomal RNA (RNR1), ribosomal RNA": "mt_rRNA_12S",
    "ref|NR_137295.1| Homo sapiens mitochondrially encoded 16S ribosomal RNA (RNR2), ribosomal RNA": "mt_rRNA_16S",
    "ref|NR_145820.1| Homo sapiens RNA, 18S ribosomal N1 (RNA18SN1), ribosomal RNA": "rRNA_18S",
}

rRNA_subtype_length = {
    "5S":121, "5.8S":159, "12S":954, "16S":1559, "18S":1869, "28S":5070, "45S":13357, "other":675
}

# Parse the unique.fa file. (Read reads ID, count, length of sequence)
def parse_unique_fa(file_path):
    """
    Parses a unique FASTA file to get read counts and lengths.
    Assumes header format: >read_id count
    """
    read_counts = {}
    read_lengths = {}
    read_sequences = {} # New: store sequences
    current_read_id = None
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                parts = line.split('\t') 
                current_read_id = parts[0][1:]
                count = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 1
                read_counts[current_read_id] = count # Record the ID and occurrence times of the sequence.
            else:
                if current_read_id:
                    read_lengths[current_read_id] = len(line) 
                    read_sequences[current_read_id] = line 
                    current_read_id = None 
    return read_counts, read_lengths, read_sequences

def get_summary_class_from_locus_and_type(locus_name, ncRNA_type):
    """
    Determines the summary class name based on the ncRNA_type and locus_name.
    """
    if ncRNA_type in ['5S', '5.8S', '12S', '16S', '18S', '28S', '45S']:
        return "rRNA_Match_Genome"
    elif ncRNA_type == "nuclear_tRNA":
        return "nuclear_tRNA_Match_Genome"
    elif ncRNA_type == "mt_RNA":
        return "mt_tRNA_Match_Genome"
    elif ncRNA_type == "yRNA":
        return "yRNA_Match_Genome"
    elif ncRNA_type == "miRNA":
        return "miRNA_Match_Genome"
    elif ncRNA_type == "snoRNA":
        return "snoRNA_Match_Genome"
    elif ncRNA_type == "snRNA":
        return "snRNA_Match_Genome" 
    elif ncRNA_type == "vaultRNA":
        return "vaultRNA_Match_Genome" 
    elif ncRNA_type == "piRNA":
        return "piRNA_Match_Genome"
    return f"Unannotated_Match_Genome" # General fallback


# Parse a separate bowtie comparison detail file.
def parse_detail_match_file(file_path, ncRNA_type, read_lengths_from_unique_fa, read_sequences_from_unique_fa):
    """
    Parses a single bowtie detail match file and extracts annotations as in annotation.pl.
    Assumes bowtie default output format: ReadID \t Counts \t Strand \t RefName \t Pos \t Sequence \t Quality \t Mismatch
    """
    read_locus_map_subset = defaultdict(set) # {read_id: set(locus_name)}
    read_to_raw_annotations_subset = defaultdict(set) # {read_id: set(raw_annotation_string)}
    read_rRNA_details_subset = defaultdict(list) # NEW: {read_id: [{'locus_name': ..., 'pos': ..., 'read_len': ..., 'strand': ..., 'full_annotation_string': ..., 'sequence': ...}, ...]}

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            parts = line.split('\t')
            if len(parts) < 6: # Need at least ReadID, Counts, Strand, RefName, Pos, Sequence
                continue
            
            read_id = parts[0]
            read_count_in_file = int(parts[1])
            strand = parts[2] 
            ref_name = parts[3] 
            pos = int(parts[4])
            seq_from_detail_file = parts[5] 
            seq = read_sequences_from_unique_fa.get(read_id, seq_from_detail_file) 
            annotation_source_field = ref_name 
            current_annotation_string = "Unannotated" 

            if ncRNA_type == "miRNA":
                # For miRNA, annotation.pl uses the first word of the RefName, e.g., "hsa-miR-1-3p"
                # Here, we use ref_name (parts[3]) as the source for annotation string.
                current_annotation_string = annotation_source_field.split()[0]

            elif ncRNA_type in ['5S', '5.8S', '12S', '16S', '18S', '28S', '45S']:
                # annotation.pl has complex regex for rRNA from its $anno field (parts[3] of their input)
                # We use the ref_name (parts[3]) as the source for annotation.
                # It distinguishes by 28S, 5.8S, 5S etc. directly from the original $anno and adds '-rRNA'.
                match_s = re.search(r'([0-9]+\.[0-9]+S|[0-9]+S)', annotation_source_field)

                rRNA_len = rRNA_subtype_length.get(ncRNA_type, 0) # Get reference length of this rRNA subtype

                # Add the classification of 5', 3', i'
                if match_s:
                    if pos == 0:
                        current_annotation_string = f"rRNA_{match_s.group(1)}_5_end"
                    elif rRNA_len > 0 and (rRNA_len == (len(seq) + pos)): 
                        current_annotation_string = f"rRNA_{match_s.group(1)}_3_end"
                    else:
                        current_annotation_string = f"rRNA_{match_s.group(1)}_inner"
                else:
                    current_annotation_string = f"rRNA_other" 
                
                # Store detailed rRNA mapping for rRNA_output.txt
                read_rRNA_details_subset[read_id].append({
                    'locus_name': ref_name,
                    'pos': pos, 
                    'read_len': len(seq), 
                    'strand': strand,
                    'full_annotation_string': current_annotation_string, 
                    'sequence': seq 
                })

            # YRNA
            elif ncRNA_type in ['yRNA']:
                match_s = re.search(r'RNY([1-5]+)', annotation_source_field)
                if match_s:
                     current_annotation_string = f"yRNA_{match_s.group(0)}"
                else:
                    current_annotation_string = "yRNA_other"

            elif ncRNA_type in ['nuclear_tRNA', 'mt_RNA']:
                # Annotation for tRNA is derived from ref_name (parts[2]), pos (parts[3]), seq (parts[4])
                ref_tRNA_len = 0
                tRNA_name = "Undet"
                tRNA_codon = "Undet"
                
                match_bp = re.search(r'([0-9]+)\s+bp', annotation_source_field)
                if match_bp:
                    ref_tRNA_len = int(match_bp.group(1))
                match_name_codon = re.search(r'([A-Za-z]+)\s+\(([A-Za-z]+)\)', annotation_source_field)
                if match_name_codon:
                    tRNA_name = match_name_codon.group(1) # Ala
                    tRNA_codon = match_name_codon.group(2) # AGC
                
                if "Homo_sapiens_tRNA" in annotation_source_field or "Homo_sapiens_mt_tRNA" in annotation_source_field:
                    prefix = "nuclear_tRNA" if ncRNA_type == 'nuclear_tRNA' else "mt_tRNA"
                    if pos == 0:
                        current_annotation_string = f"{prefix}-{tRNA_name}-{tRNA_codon}_5_end"
                    elif (pos + len(seq) == ref_tRNA_len) and ref_tRNA_len > 0 : 
                         current_annotation_string = f"{prefix}-{tRNA_name}-{tRNA_codon}_3_end"
                    elif (pos + len(seq) + 3 == ref_tRNA_len) and ref_tRNA_len > 0: # +3 for CCA at 3' end
                        current_annotation_string = f"{prefix}-{tRNA_name}-{tRNA_codon}_3_end_CCA"
                    else:
                        current_annotation_string = f"{prefix}-{tRNA_name}-{tRNA_codon}_inner"
                else:
                    current_annotation_string = f"{tRNA_name}-{tRNA_codon}" # Fallback
            
            elif ncRNA_type == 'piRNA':
                current_annotation_string = "piRNA"
            
            elif ncRNA_type == 'vaultRNA':
                match_s = re.search(r'VTRNA([0-9]+)-([0-9]+)', annotation_source_field)
                if match_s:
                    current_annotation_string = f"VTRNA{match_s.group(1)}-{match_s.group(2)}"
                else: 
                    current_annotation_string = "other_vtRNA"
            
            elif ncRNA_type == 'snoRNA':
                current_annotation_string = "snoRNA"
            
            elif ncRNA_type == 'snRNA':
                current_annotation_string = "snRNA"

            # Fallback for other ncRNA types not explicitly handled above (e.g., snoRNA, snRNA, vaultRNA, yRNA without specific regex)
            else: 
                current_annotation_string = ncRNA_type

            read_locus_map_subset[read_id].add(ref_name) # Locus for weighted allocation
            read_to_raw_annotations_subset[read_id].add(current_annotation_string) # For output.txt

    return read_locus_map_subset, read_to_raw_annotations_subset, read_rRNA_details_subset

def collect_all_mappings(sample_name, output_dir, read_lengths_from_unique_fa, read_sequences_from_unique_fa):
    """
    Collects all mapping information from various detail match files.
    Returns:
        all_read_locus_map: {read_id: set(locus_name)}
        all_locus_to_class_map: {locus_name: summary_class_string}
        all_read_to_annotations: {read_id: set(raw_annotation_string)}
        all_rRNA_detailed_mappings: {read_id: list of rRNA mapping dicts}
    """
    all_read_locus_map = defaultdict(set)
    all_locus_to_class_map = {}
    all_read_to_annotations = defaultdict(set)
    all_rRNA_detailed_mappings = defaultdict(list) # NEW!

    # List of tuples: (file_suffix, ncRNA_type_for_parsing)
    ncRNA_types_to_process = [
        ('miRNA_ncRNA', 'miRNA'), 
        ('5S_rRNA', '5S'), ('5.8S_rRNA', '5.8S'), ('12S_rRNA', '12S'), ('16S_rRNA', '16S'), ('18S_rRNA', '18S'), ('28S_rRNA', '28S'), ('45S_rRNA', '45S'),
        ('nuclear_tRNA_ncRNA', 'nuclear_tRNA'), ('mt_RNA_ncRNA', 'mt_RNA'), 
        ('yRNA_ncRNA', 'yRNA'), 
        ('snoRNA_ncRNA', 'snoRNA'), 
        ('snRNA_ncRNA', 'snRNA'), 
        ('vaultRNA_ncRNA', 'vaultRNA'), 
        ('piRNA_ncRNA', 'piRNA')
    ]

    for file_info in ncRNA_types_to_process:
        file_suffix, ncRNA_type = file_info
        detail_file_path = os.path.join(output_dir, "mapped", f"{sample_name}_detail_match_{file_suffix}")
        
        if os.path.exists(detail_file_path):
            print(f"Parsing {detail_file_path} for {ncRNA_type}...")
            read_locus_map_subset, read_to_raw_annotations_subset, read_rRNA_details_subset = \
                parse_detail_match_file(detail_file_path, ncRNA_type, read_lengths_from_unique_fa, read_sequences_from_unique_fa) # Pass read_sequences_from_unique_fa

            for read_id, loci in read_locus_map_subset.items():
                all_read_locus_map[read_id].update(loci)
                # Populate all_locus_to_class_map here
                for locus_name in loci:
                    all_locus_to_class_map[locus_name] = get_summary_class_from_locus_and_type(locus_name, ncRNA_type)
            
            for read_id, annotations_set in read_to_raw_annotations_subset.items():
                all_read_to_annotations[read_id].update(annotations_set)

            # Aggregate rRNA details (only if it's an rRNA type)
            if ncRNA_type in ['5S', '5.8S', '12S', '16S', '18S', '28S', '45S']:
                for read_id, details_list in read_rRNA_details_subset.items():
                    all_rRNA_detailed_mappings[read_id].extend(details_list) # Use extend for list of dicts

    return all_read_locus_map, all_locus_to_class_map, all_read_to_annotations, all_rRNA_detailed_mappings

# Function to realize weighted distribution
def calculate_weighted_counts(read_counts, read_lengths_param, all_read_locus_map): # NEW PARAS
    """
    Calculates weighted counts for each locus and their length distribution.
    Arguments:
        read_counts (dict): {read_id: original_count}
        read_lengths_param (dict): {read_id: length} (now passed as parameter)
        all_read_locus_map (defaultdict(set)): {read_id: set(locus_name)}
    Returns:
        final_locus_counts (defaultdict(float)): {locus: allocated_count}
        final_locus_length_dist (defaultdict(defaultdict(float))): {locus: {length: allocated_count}}
    """
    final_locus_counts = defaultdict(float)
    final_locus_length_dist = defaultdict(lambda: defaultdict(float))
    
    # Step 1: Calculate initial unique abundances and populate final counts for unique reads
    unique_locus_abundance = defaultdict(float) # Using float for consistency with later sums

    for read_id, loci in all_read_locus_map.items():
        if len(loci) == 1:
            locus = list(loci)[0]
            current_read_count = float(read_counts.get(read_id, 1))
            unique_locus_abundance[locus] += current_read_count
            final_locus_counts[locus] += current_read_count
            # For length distribution, directly add to the unique locus's length distribution
            final_locus_length_dist[locus][read_lengths_param.get(read_id, 0)] += current_read_count
        else:
            # Ensure multi-mapped loci are initialized in final_locus_counts and length_dist
            for locus in loci:
                final_locus_counts.setdefault(locus, 0.0)
                final_locus_length_dist.setdefault(locus, defaultdict(float))


    # Step 2: Distribute multi-mapping reads using weighted approach
    for read_id, loci in all_read_locus_map.items():
        if len(loci) > 1: # Only process multi-mapping reads
            current_read_count = float(read_counts.get(read_id, 1))
            read_len = read_lengths_param.get(read_id, 0) # Use parameter
            
            # Calculate total unique abundance for all loci this read maps to
            total_unique_abundance_of_mapped_loci = sum(unique_locus_abundance.get(locus, 0.0) for locus in loci)

            if total_unique_abundance_of_mapped_loci == 0:
                # If all loci have zero unique mapping reads, distribute evenly
                share = current_read_count / len(loci)
                for locus in loci:
                    final_locus_counts[locus] += share
                    final_locus_length_dist[locus][read_len] += share
            else:
                # Distribute based on unique abundance as weights
                for locus in loci:
                    weight = unique_locus_abundance.get(locus, 0.0) / total_unique_abundance_of_mapped_loci
                    allocated_count = current_read_count * weight
                    final_locus_counts[locus] += allocated_count
                    final_locus_length_dist[locus][read_len] += allocated_count
    
    return final_locus_counts, final_locus_length_dist

# Calculate Unannotated reads 
def parse_unannotated_fa(file_path):
    """
    Parses the unannotated FASTA file for total count and length distribution.
    Assumes header format: >read_id\tcount
    """
    unannotated_total_count = 0.0
    unannotated_length_distribution = defaultdict(float) # {length: count}
    
    if not os.path.exists(file_path):
        print(f"Warning: Unannotated FASTA file not found: {file_path}", file=sys.stderr)
        return unannotated_total_count, unannotated_length_distribution

    with open(file_path, 'r') as f:
        current_read_id = None # Need to track count per read_id, as sequences follow header
        current_read_count = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                parts = line.split('\t') # Assuming tab-separated header
                current_read_id = parts[0][1:]
                current_read_count = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 1
                unannotated_total_count += current_read_count
            else:
                if current_read_id: # Only process if a header was just read
                    seq_length = len(line)
                    unannotated_length_distribution[seq_length] += current_read_count # Use the 'count' from the header
                    current_read_id = None # Reset for next read
            
    return unannotated_total_count, unannotated_length_distribution

# Output: output.txt
def write_output_txt(file_path, read_counts, read_lengths, read_sequences, read_match_genome_status, all_read_to_annotations): # NEW ARGS
    """
    Writes the output.txt file based on annotation.pl's format.
    ID\tSequence\tLength\tReads\tMatch_Genome\tAnnotation
    """
    with open(file_path, "w") as f:
        f.write("ID\tSequence\tLength\tReads\tMatch_Genome\tAnnotation\n")
        
        # Iterate through all read_ids from the original unique.fa
        for read_id in sorted(read_counts.keys()):
            sequence = read_sequences.get(read_id, "")
            length = read_lengths.get(read_id, 0)
            reads_count = read_counts.get(read_id, 0)
            match_genome = read_match_genome_status.get(read_id, "NO")
            
            # Combine annotations, if any
            annotations = all_read_to_annotations.get(read_id, set())
            
            if annotations:
                annotation_str = ";".join(sorted(list(annotations))) # Sort for consistent output
            else:
                annotation_str = "NO_Annotation"
            
            # If it didn't match genome, force annotation to "NO_Annotation"
            if match_genome == "NO":
                annotation_str = "NO_Annotation"

            f.write(f"{read_id}\t{sequence}\t{length}\t{reads_count}\t{match_genome}\t{annotation_str}\n")

# Output: summary.txt
def write_summary_txt(file_path, total_clean_reads, final_locus_counts, all_locus_to_class_map, unannotated_total_count, locus_name_map):
    """
    Writes the summary.txt file (aggregated counts by class).
    """
    summary_class_counts = defaultdict(float)
    total_annotated_reads = 0.0
    
    for locus, count in final_locus_counts.items():
        summary_class = all_locus_to_class_map.get(locus, "Unknown_Class_Match_Genome")
        summary_class_counts[summary_class] += count
        total_annotated_reads += count

    total_mapped_reads = total_annotated_reads + unannotated_total_count

    with open(file_path, "w") as f:
        f.write("Class\tSub_Class\tReads\n")
        f.write(f"Clean_Reads\t-\t{total_clean_reads}\n")
        f.write(f"Match_Genome\t-\t{total_mapped_reads:.2f}\n")
        
        for class_name in sorted(summary_class_counts.keys()):
            f.write(f"{class_name}\t-\t{summary_class_counts[class_name]:.2f}\n")
            for locus, count in sorted(final_locus_counts.items()):
                if class_name == "rRNA_Match_Genome" and all_locus_to_class_map.get(locus) == class_name:
                    simplified_locus = locus_name_map.get(locus, locus) 
                    f.write(f"{class_name}\t{simplified_locus}\t{count:.2f}\n")
        
        f.write(f"Unannotated_Match_Genome\t-\t{unannotated_total_count:.2f}\n")

# Output: length_distribution.txt
def write_length_distribution_txt(file_path, read_lengths, read_counts, all_read_locus_map, final_locus_length_dist, all_locus_to_class_map, unannotated_length_distribution):
    """
    Writes the length_distribution.txt file.
    """
    class_length_dist = defaultdict(lambda: defaultdict(float))

    for read_id in read_lengths:
        length = read_lengths[read_id]
        count = read_counts.get(read_id, 1)
        class_length_dist["Clean_Reads"][length] += count

    for length, count in unannotated_length_distribution.items():
        class_length_dist["Unannotated_Match_Genome"][length] += count 

    for locus, length_dist_map in final_locus_length_dist.items():
        summary_class = all_locus_to_class_map.get(locus, "Unknown_Class_Match_Genome")
        for length, count in length_dist_map.items():
            class_length_dist[summary_class][length] += count

    with open(file_path, "w") as f:
        f.write("Class\tLength\tReads\n")
        for class_name in sorted(class_length_dist.keys()):
            for length in sorted(class_length_dist[class_name].keys()):
                count = class_length_dist[class_name][length]
                f.write(f"{class_name}\t{length}\t{count:.2f}\n")

# Output: rRNA_output.txt
def write_rRNA_output_txt(output_file, read_counts_dict, all_rRNA_detailed_mappings):
    """
    Generates the rRNA_output.txt file with weighted counts for multi-mapping rRNAs.
    """
    with open(output_file, 'w') as out_f:
        out_f.write("Sequence\tOrigin\tSubtype\tFragment\tLength\tStart\tEnd\tCounts\n")

        unique_rRNA_locus_segment_counts = defaultdict(float)

        for read_id, mappings in all_rRNA_detailed_mappings.items():
            original_count = read_counts_dict.get(read_id, 0) # collapse_read counts
            if original_count == 0:
                continue

            unique_mapping_instances_for_read = set()
            for mapping in mappings:
                unique_mapping_instances_for_read.add((mapping['locus_name'], mapping['pos'], mapping['read_len']))
            
            if len(unique_mapping_instances_for_read) == 1:
                mapping = mappings[0]
                locus_key = (mapping['locus_name'], mapping['pos'], mapping['read_len'])
                unique_rRNA_locus_segment_counts[locus_key] += original_count

        processed_output_lines = []
        for read_id, mappings in all_rRNA_detailed_mappings.items():
            original_read_count = read_counts_dict.get(read_id, 0)
            if original_read_count == 0:
                continue # Skip reads with no count

            grouped_mappings = defaultdict(list)
            for mapping in mappings:
                locus_key = (mapping['locus_name'], mapping['pos'], mapping['read_len'])
                grouped_mappings[locus_key].append(mapping) 

            if len(grouped_mappings) == 1:
                locus_key = list(grouped_mappings.keys())[0]
                representative_mapping = grouped_mappings[locus_key][0] 
                allocated_count = original_read_count
                processed_output_lines.append(
                    format_rRNA_output_line(representative_mapping, allocated_count)
                )
            else:
                total_unique_abundance_for_multi_mapper = 0
                for locus_key in grouped_mappings.keys():
                    total_unique_abundance_for_multi_mapper += unique_rRNA_locus_segment_counts[locus_key]
                
                if total_unique_abundance_for_multi_mapper > 0:
                    # Weighted allocation based on unique counts of the loci involved
                    for locus_key, specific_mappings in grouped_mappings.items():
                        representative_mapping = specific_mappings[0] 
                        locus_unique_abundance = unique_rRNA_locus_segment_counts[locus_key]
                        
                        allocated_count = original_read_count * (locus_unique_abundance / total_unique_abundance_for_multi_mapper)
                        if allocated_count > 0: # Only add if allocated count is greater than 0
                            processed_output_lines.append(
                                format_rRNA_output_line(representative_mapping, allocated_count)
                            )
                else:
                    # Fallback to even distribution if no unique abundance (all loci have 0 unique reads)
                    # This happens if all reads mapping to these loci are themselves multi-mappers
                    allocated_count_per_locus = original_read_count / len(grouped_mappings)
                    for locus_key, specific_mappings in grouped_mappings.items():
                        representative_mapping = specific_mappings[0]
                        if allocated_count_per_locus > 0:
                            processed_output_lines.append(
                                format_rRNA_output_line(representative_mapping, allocated_count_per_locus)
                            )
        
        # Sort lines for consistent output. Sort by Sequence, then Start, then End for deterministic order.
        processed_output_lines.sort(key=lambda x: (x.split('\t')[0], int(x.split('\t')[5]), int(x.split('\t')[6])))
        for line in processed_output_lines:
            out_f.write(line + '\n')

def format_rRNA_output_line(mapping_detail, allocated_count):
    """Helper function to format a single line for rRNA_output.txt"""
    seq = mapping_detail['sequence']
    locus_name = mapping_detail['locus_name']
    pos = mapping_detail['pos']
    read_len = mapping_detail['read_len']
    full_annotation_string = mapping_detail['full_annotation_string'] # e.g., "rRNA_28S_5_end"

    # Determine Origin
    if '12S' in locus_name or '16S' in locus_name or 'mt_rRNA' in full_annotation_string: 
        origin = "mito"
    else:
        origin = "nuclear-encoded"

    # Determine Subtype (e.g., 5S, 12S, 28S) from full_annotation_string
    subtype_match = re.search(r'rRNA_([0-9]+\.[0-9]+S|[0-9]+S|other)', full_annotation_string)
    subtype = subtype_match.group(1) if subtype_match else "Unknown"
    if "mt_rRNA_12S" in full_annotation_string: subtype = "12S"
    if "mt_rRNA_16S" in full_annotation_string: subtype = "16S"

    # Determine Fragment (e.g., 5_end, inner, 3_end)
    fragment_match = re.search(r'(5_end|3_end|inner)', full_annotation_string) # Include CCA_end from tRNA logic if it sneaks in
    fragment = fragment_match.group(1) if fragment_match else "Full" # Fallback to Full or Other

    start = pos + 1 
    end = pos + read_len 

    return f"{seq}\t{origin}\t{subtype}\t{fragment}\t{read_len}\t{start}\t{end}\t{allocated_count:.4f}"


def write_rRNA_length_distribution_txt(file_path, final_read_counts_per_locus_type, read_lengths_from_unique_fa, all_rRNA_detailed_mappings):
    """
    Writes the length distribution of rRNA subtypes to a file.

    Args:
        file_path (str): The path to the output file (e.g., sample_rRNA_length_distribution.txt).
        final_read_counts_per_locus_type (dict): Dictionary mapping read_id to a dictionary of (locus_name, ncRNA_type) to allocated_count.
        read_lengths_from_unique_fa (dict): Dictionary mapping read_id to its length.
        all_rRNA_detailed_mappings (dict): Dictionary mapping read_id to a list of detailed rRNA mappings.
    """
    rRNA_subtype_length_dist = {
        "5S": {}, "5.8S": {}, "12S": {}, "16S": {}, "18S": {}, "28S": {}, "45S": {}
    }

    for read_id, locus_type_counts in final_read_counts_per_locus_type.items():
        read_length = read_lengths_from_unique_fa.get(read_id)
        if read_length is None:
            continue 

        for (locus_name, ncRNA_type), allocated_count in locus_type_counts.items():
            if ncRNA_type == "rRNA":
                matching_rRNA_detail = None
                for detail in all_rRNA_detailed_mappings.get(read_id, []):
                    if detail['locus_name'] == locus_name:
                        matching_rRNA_detail = detail
                        break
                
                if matching_rRNA_detail:
                    fragment = matching_rRNA_detail['fragment'] 
                    subtype = fragment.split('_rRNA')[0]
                    
                    if subtype in rRNA_subtype_length_dist:
                        rRNA_subtype_length_dist[subtype][read_length] = \
                            rRNA_subtype_length_dist[subtype].get(read_length, 0) + allocated_count

    # Write to file
    with open(file_path, 'w') as f:
        f.write("rRNA_subtype\tLength\tCounts\n")
        
        sorted_subtypes = sorted(rRNA_subtype_length_dist.keys())
        for subtype in sorted_subtypes:
            sorted_lengths = sorted(rRNA_subtype_length_dist[subtype].keys())
            for length in sorted_lengths:
                count = rRNA_subtype_length_dist[subtype][length]
                f.write(f"{subtype}\t{length}\t{round(count, 3)}\n")


def main(sample_name, output_dir, ncrna_flag):
    """
    Main function to orchestrate the RNA-seq processing.
    """
    # Input file paths
    unique_fa_file = os.path.join(output_dir, "mapped", f"{sample_name}_unique.fa")
    
    if ncrna_flag == "true": 
        unannotated_piRNA_fa_file = os.path.join(output_dir, "mapped",f"{sample_name}_MG_unmatch_piRNA_ncRNA.fa")
    else:
        unannotated_piRNA_fa_file = os.path.join(output_dir, "mapped",f"{sample_name}_MG_unmatch_45S_rRNA.fa")


    # Parse unique FASTA for read counts, lengths, and sequences
    read_counts, read_lengths, read_sequences = parse_unique_fa(unique_fa_file) # NEW!
    total_clean_reads = sum(read_counts.values())

    # Collect all annotated mappings and raw annotations
    all_read_locus_map, all_locus_to_class_map, all_read_to_annotations, all_rRNA_detailed_mappings = \
        collect_all_mappings(sample_name, output_dir, read_lengths, read_sequences) # NEW!
    
    # Determine Match_Genome status for all reads
    read_match_genome_status = {}
    for read_id in read_counts: # Iterate all reads from unique.fa
        if read_id in all_read_locus_map: # If it has any mapped locus
            read_match_genome_status[read_id] = "Yes"
        else:
            read_match_genome_status[read_id] = "NO"

    # Calculate weighted counts and locus-level length distributions
    final_locus_counts, final_locus_length_dist = calculate_weighted_counts(read_counts, read_lengths, all_read_locus_map) # NEW PARAM!

    # Parse unannotated FASTA
    unannotated_total_count, unannotated_length_distribution = parse_unannotated_fa(unannotated_piRNA_fa_file)

    # Output file paths
    output_txt_file = os.path.join(output_dir, "stats", f"{sample_name}_output.txt")
    summary_txt_file = os.path.join(output_dir, "stats", f"{sample_name}_summary.txt")
    length_distribution_txt_file = os.path.join(output_dir, "stats", f"{sample_name}_length_distribution.txt")

    # Write output files
    write_output_txt(output_txt_file, read_counts, read_lengths, read_sequences, read_match_genome_status, all_read_to_annotations)
    write_summary_txt(summary_txt_file, total_clean_reads, final_locus_counts, all_locus_to_class_map, unannotated_total_count, locus_name_map)
    write_length_distribution_txt(length_distribution_txt_file, read_lengths, read_counts, all_read_locus_map, final_locus_length_dist, all_locus_to_class_map, unannotated_length_distribution)

    # Generate rRNA_output.txt (NEW!)
    rRNA_output_file = os.path.join(output_dir, "stats", f"{sample_name}_rRNA_output.txt")
    print(f"Generating {rRNA_output_file}...")
    write_rRNA_output_txt(rRNA_output_file, read_counts, all_rRNA_detailed_mappings)
    print(f"Successfully generated {rRNA_output_file}.")

    print(f"Processing for sample {sample_name} completed. Results in {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_rna_seq.py <sample_name> <output_directory> <ncrna_flag>")
        sys.exit(1)
    
    sample_name = sys.argv[1]
    output_dir = sys.argv[2]
    ncrna_flag = sys.argv[3]

    main(sample_name, output_dir, ncrna_flag) 
