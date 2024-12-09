import pandas as pd
import primer3
import copy
import csv
from itertools import combinations
from collections import Counter
import numpy as np
import difflib as dl
import matplotlib.pyplot as plt
import statistics

def CAR_sense_strand_creation():
    copied_whole_seq = "ATGGCCCTGCCGGTGACGGCCCTGCTGCTGCCCCTGGCGCTGCTTTTGCATGCGGCTAGGCCTGAGATCGTGATGACCCAGAGTCCCGCCACCCTTTCCGTAAGCCCCGGGGAGCGGGCTACACTGTCATGCCGAGCGTCACAGAGCGTCGACAACAACCTGGTGTGGTACCAGCAGAAGCCTGGACAGGCCCCGCGCCTGCTCATTTACGGCGCTTCTACCCGTGCTACCGGCATCCCAGCTCGTTTTTCAGGCTCCGGATCCGGGACCGAGTTCACCCTGACGATCTCTTCTTTGCAGAGCGAAGATTTTGCGGTTTATTACTGTCAACAGTACAACAATTGGCCGCCGTGGACCTTCGGTCAGGGCACCAAGGTGGAGATCAAGCGCGGGGGCGGTGGCTCCGGGGGTGGCGGGAGCGGAGGCGGTGGTAGCCAGGTGCAGCTGCAGGAGAGCGGCCCTGGCCTCGTGAAGCCCTCTCAGACTCTTTCGCTAACCTGCACTGTCTCGGGCGGGTCCATCAGCTCCTTCAACTACTACTGGAGTTGGATCCGCCACCACCCCGGCAAAGGCCTGGAGTGGATTGGCTACATCTACTACTCTGGTTCCACCTACTCCAATCCATCTCTGAAATCCCGCGTGACTATTTCTGTCGATACCTCCAAGAACCAGTTCTCCCTGACCCTGAGTTCCGTGACTGCTGCGGACACCGCCGTGTACTACTGTGCCCGCGGCTACAACTGGAACTATTTCGACTATTGGGGCCAGGGAACTCTCGTCACCGTGTCGAGCGCAGAACAGAAGCTCATCTCGGAGGAGGACCTGGCCGCAGCCACCACAACCccagcgccgcgaccaccaacaccggcgcccaccatcgcgtcgcagcccctgtccctgcgcccagaggcgtgccggccagcggcggggggcgcagtgcacacgagggggctggacttcgcctgtgatatctacatctgggcgcccttggccgggacttgtggggtccttctcctgtcactggttatcaccctttactgcaaacggggcagaaagaaactcctgtatatattcaaacaaccatttatgagaccagtacaaactactcaagaggaagatggctgtagctgccgatttccagaagaagaagaaggaggatgtgaactgagagtgaagttcagcaggagcgcagacgcccccgcgtaccagcagggccagaaccagctctataacgagctcaatctaggacgaagagaggagtacgatgttttggacaagagacgtggccgggaccctgagatggggggaaagccgagaaggaagaaccctcaggaaggcctgtacaatgaactgcagaaagataagatggcggaggcctacagtgagattgggatgaaaggcgagcgccggaggggcaaggggcacgatggcctttaccagggtctcagtacagccaccaaggacacctacgacgcccttcacatgcaggccctgccccctcgc"

    cd8a_leader = "ATGGCCCTGCCGGTGACGGCCCTGCTGCTGCCCCTGGCGCTGCTTTTGCATGCGGCTAGGCCT"
    vi_gpnmb = "GAGATCGTGATGACCCAGAGTCCCGCCACCCTTTCCGTAAGCCCCGGGGAGCGGGCTACACTGTCATGCCGAGCGTCACAGAGCGTCGACAACAACCTGGTGTGGTACCAGCAGAAGCCTGGACAGGCCCCGCGCCTGCTCATTTACGGCGCTTCTACCCGTGCTACCGGCATCCCAGCTCGTTTTTCAGGCTCCGGATCCGGGACCGAGTTCACCCTGACGATCTCTTCTTTGCAGAGCGAAGATTTTGCGGTTTATTACTGTCAACAGTACAACAATTGGCCGCCGTGGACCTTCGGTCAGGGCACCAAGGTGGAGATCAAGCGC"
    g4s_linker = "GGGGGCGGTGGCTCCGGGGGTGGCGGGAGCGGAGGCGGTGGTAGC"
    vh_gpnmb = "CAGGTGCAGCTGCAGGAGAGCGGCCCTGGCCTCGTGAAGCCCTCTCAGACTCTTTCGCTAACCTGCACTGTCTCGGGCGGGTCCATCAGCTCCTTCAACTACTACTGGAGTTGGATCCGCCACCACCCCGGCAAAGGCCTGGAGTGGATTGGCTACATCTACTACTCTGGTTCCACCTACTCCAATCCATCTCTGAAATCCCGCGTGACTATTTCTGTCGATACCTCCAAGAACCAGTTCTCCCTGACCCTGAGTTCCGTGACTGCTGCGGACACCGCCGTGTACTACTGTGCCCGCGGCTACAACTGGAACTATTTCGACTATTGGGGCCAGGGAACTCTCGTCACCGTGTCGAGCGCA"
    myc = "GAACAGAAGCTCATCTCGGAGGAGGACCTG"
    post_linker = "GCCGCAGCCACCACAACC"
    cd8a = "ccagcgccgcgaccaccaacaccggcgcccaccatcgcgtcgcagcccctgtccctgcgcccagaggcgtgccggccagcggcggggggcgcagtgcacacgagggggctggacttcgcctgtgatatctacatctgggcgcccttggccgggacttgtggggtccttctcctgtcactggttatcaccctttactgc"
    isd_41bb = "aaacggggcagaaagaaactcctgtatatattcaaacaaccatttatgagaccagtacaaactactcaagaggaagatggctgtagctgccgatttccagaagaagaagaaggaggatgtgaa"
    cd3z = "ctgagagtgaagttcagcaggagcgcagacgcccccgcgtaccagcagggccagaaccagctctataacgagctcaatctaggacgaagagaggagtacgatgttttggacaagagacgtggccgggaccctgagatggggggaaagccgagaaggaagaaccctcaggaaggcctgtacaatgaactgcagaaagataagatggcggaggcctacagtgagattgggatgaaaggcgagcgccggaggggcaaggggcacgatggcctttaccagggtctcagtacagccaccaaggacacctacgacgcccttcacatgcaggccctgccccctcgc"

    # Construct the whole sequence, with the N's as the domain junctions.
    constructed_whole_seq = cd8a_leader + "N" + vi_gpnmb + "N" + g4s_linker + "N" + vh_gpnmb + "N" + myc + "N" + post_linker + "N" + cd8a + "N" + isd_41bb + "N" + cd3z

    constructed_whole_seq = constructed_whole_seq.lower()
    constructed_whole_seq = constructed_whole_seq.replace("n", "N")

    return constructed_whole_seq

def strand_creation(sense_strand):
    replacements = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'N': 'N'}
    replaced_chars = [replacements.get(char, char) for char in sense_strand]
    antisense_stand = ''.join(replaced_chars)

    return antisense_stand

def probes_from_sequence(reverse_antisense_strand, total_probe_length):
    possible_probes = []

    for i in range(len(reverse_antisense_strand)-total_probe_length+1):
        N_count = reverse_antisense_strand[i:i+total_probe_length].count("N")

        # If there's an "N" in the probe, we add N number of extra nucleotide to the probe. The N is indicative of the probe spanning a domain junction.
        if N_count >= 1:
            probe_seq = reverse_antisense_strand[i:i+total_probe_length+N_count]
        else:
            probe_seq = reverse_antisense_strand[i:i+total_probe_length]
        possible_probes.append(probe_seq)

    return possible_probes

def split_probes(probes):
    probe_table = dict(zip(range(0, len(probes)), [probes]))
    probe_table = pd.DataFrame(probe_table)
    probe_table.columns= ["probe_sequence"]

    LHS_probe_list = []
    RHS_probe_list = []
    domain_spanning_list = []

    for row, index in probe_table.iterrows():
        temp = index["probe_sequence"]

        if "N" in index["probe_sequence"]:
            domain_spanning_list.append( temp.count("N") )
            temp = temp.replace("N", "")

            LHS_probe_list.append(temp[:25])
            RHS_probe_list.append(temp[25:])
        elif "N" not in index["probe_sequence"]:
            LHS_probe_list.append(temp[:25])
            RHS_probe_list.append(temp[25:])
            domain_spanning_list.append(0)

        probe_table.loc[row]["domain_spanning"] = 1 if "N" in index["probe_sequence"] else 0

    probe_table["LHS_probe"] = LHS_probe_list
    probe_table["RHS_probe"] = RHS_probe_list
    probe_table["domain_spanning"] = domain_spanning_list
    
    return probe_table

def longest_stretch(string):
    max_stretch = 0
    current_stretch = 1
    
    for i in range(1, len(string)):
        if string[i] == string[i-1]:
            current_stretch += 1
        else:
            max_stretch = max(max_stretch, current_stretch)
            current_stretch = 1
            
    max_stretch = max(max_stretch, current_stretch)
    return max_stretch

def calculate_probe_values(probes_table):
    # Get the GC content of the probes.
    probes_table["LHS_GC"] = probes_table["LHS_probe"].apply(lambda x: (((x.count("g") + x.count("c"))/25)*100) )
    probes_table["RHS_GC"] = probes_table["RHS_probe"].apply(lambda x: (((x.count("g") + x.count("c"))/25)*100) )
    probes_table["probe_GC_diff"] = probes_table["LHS_GC"] - probes_table["RHS_GC"]
    probes_table["probe_GC"] = probes_table["LHS_GC"] + probes_table["RHS_GC"]

    # Get the longest strech of consecutive nucleotides in the probes.
    probes_table["LHS_homopolymer"] = probes_table["LHS_probe"].apply(longest_stretch)
    probes_table["RHS_homopolymer"] = probes_table["RHS_probe"].apply(longest_stretch)
    probes_table["probe_homopolymer"] = probes_table["LHS_homopolymer"] + probes_table["RHS_homopolymer"]

    probes_table["LHS_homopolymer_capped"] = probes_table["LHS_homopolymer"]
    probes_table["RHS_homopolymer_capped"] = probes_table["RHS_homopolymer"]

    probes_table.loc[probes_table["LHS_homopolymer_capped"] <= 4, "LHS_homopolymer_capped"] = 0
    probes_table.loc[probes_table["LHS_homopolymer_capped"] >= 5, "LHS_homopolymer_capped"] = 1

    probes_table.loc[probes_table["RHS_homopolymer_capped"] <= 4, "RHS_homopolymer_capped"] = 0
    probes_table.loc[probes_table["RHS_homopolymer_capped"] >= 5, "RHS_homopolymer_capped"] = 1

    # Get the GC content to homopolymer ratio.
    probes_table["GC_homopolymer"] = probes_table["probe_GC"] / probes_table["probe_homopolymer"]

    # Calculate the melting temperature of the probes.
    probes_table["probe_tm"] = probes_table["probe_sequence"].apply(lambda x: primer3.calc_tm(x.upper()))
    probes_table["LHS_tm"] = probes_table["LHS_probe"].apply(lambda x: primer3.calc_tm(x.upper()))
    probes_table["RHS_tm"] = probes_table["RHS_probe"].apply(lambda x: primer3.calc_tm(x.upper()))

    # Calculate the hairpin structure of the probes.
    probes_table["probe_hairpin"] = probes_table["probe_sequence"].apply(lambda x: primer3.calc_hairpin(x.upper()).structure_found )
    probes_table["LHS_hairpin"] = probes_table["LHS_probe"].apply(lambda x: primer3.calc_hairpin(x.upper()).structure_found )
    probes_table["RHS_hairpin"] = probes_table["RHS_probe"].apply(lambda x: primer3.calc_hairpin(x.upper()).structure_found )

    probes_table["probe_hairpin_tm"] = probes_table["probe_sequence"].apply(lambda x: primer3.calc_hairpin(x.upper()).tm )
    probes_table["LHS_hairpin_tm"] = probes_table["LHS_probe"].apply(lambda x: primer3.calc_hairpin(x.upper()).tm )
    probes_table["RHS_hairpin_tm"] = probes_table["RHS_probe"].apply(lambda x: primer3.calc_hairpin(x.upper()).tm )
    
    # Calculate the number of tandem repeats in the probes.
    for index, row in probes_table.iterrows():
        LHS_probe = probes_table["LHS_probe"][index]
        RHS_probe = probes_table["RHS_probe"][index]

        LHS_probe_repeats = []
        RHS_probe_repeats = []

        for i, j in combinations(range(25 + 1), 2):
            LHS_probe_repeats.append(LHS_probe[i:j])
            RHS_probe_repeats.append(RHS_probe[i:j])
        LHS_output = [x for x in LHS_probe_repeats if len(x) >= 2]
        RHS_output = [x for x in RHS_probe_repeats if len(x) >= 2]

        LHS_repeats = Counter(LHS_output)
        RHS_repeats = Counter(RHS_output)

        LHS_repeats_subset = {k:v for (k,v) in LHS_repeats.items() if v > 1}
        RHS_repeats_subset = {k:v for (k,v) in RHS_repeats.items() if v > 1}

        probes_table.loc[index, "LHS_repeats"] = str(LHS_repeats_subset)
        probes_table.loc[index, "RHS_repeats"] = str(RHS_repeats_subset)
    
    return probes_table

def filter_probes(probes_table, filter = "default"):
    filtered_table = probes_table

    if filter == "LHS_T_end":
        filtered_table = filtered_table[probes_table["LHS_probe"].str.endswith("t")]
        print("Probes with an LHS 3' of T:", len(filtered_table))

    elif filter == "GC_content":
        filtered_table = filtered_table[(filtered_table["LHS_GC"] >= 44) & (filtered_table["RHS_GC"] >= 44) & 
                                        (filtered_table["LHS_GC"] <= 72) & (filtered_table["RHS_GC"] <= 72)]
        print("Probes with between 44-72'%' GC content:", len(filtered_table))
   
    elif filter == "homopolymer_repeat":
        filtered_table = filtered_table[filtered_table["LHS_homopolymer_capped"] == 0]
        filtered_table = filtered_table[filtered_table["RHS_homopolymer_capped"] == 0]

        print("Probes with less than or equal to 4 homopolymer repeats:", len(filtered_table))

    elif filter == "strict":
        filtered_table = filtered_table[filtered_table["LHS_homopolymer"] <= 3] 
        filtered_table = filtered_table[filtered_table["RHS_homopolymer"] <= 3]
        filtered_table = filtered_table[filtered_table["probe_GC"] >= 116]

    return filtered_table

def probe_fasta(probes_table):
    probes_fasta = probes_table.apply(lambda row: f">probe_{row.name}\n{row['probe_sequence']}", axis=1)

    if len(probes_table) == probes:
        probes_fasta.to_csv("CAR/CAR_raw_probe_fasta.fasta", index=False, header=False, sep="\n", quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ")
    else:
        probes_fasta.to_csv("CAR/CAR_filtered_probe_fasta.fasta", index=False, header=False, sep="\n", quoting=csv.QUOTE_NONE, quotechar="",  escapechar=" ")

def blast_probes(probes_table, blast_db):
    blast_results = pd.read_csv(blast_db, header = None)
    blast_results.columns = ["probe_ID", "subject_ID", "percent_identity", "length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]

    dropped = []

    for index, row in blast_results.iterrows():
        if row["mismatches"] <= 5 and row["length"] > 45:
            probe_index = int(row["probe_ID"].split("_")[1])
            if probe_index in probes_table.index:
                print("Probe %s in probes table" % probe_index)
                probes_table = probes_table.drop(index = probe_index, inplace = False)
                dropped.append(probe_index)

    return probes_table

def blast_probes_mismatches(probes_table, blast_db, mismatch_col = "probe_mismatches", probe_length_col = "probe_length"):
    blast_results = pd.read_csv(blast_db, header = None)
    blast_results.columns = ["probe_ID", "subject_ID", "percent_identity", "length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    blast_results["probe"] = blast_results["probe_ID"].apply(lambda x: x.split("_")[1])

    probes_table[mismatch_col] = ""
    probes_table[probe_length_col] = ""

    for index, row in probes_table.iterrows():
        mismatch_probes = list(blast_results[blast_results["probe"] == str(row.name)]["mismatches"])
        length_probes = list(blast_results[blast_results["probe"] == str(row.name)]["length"])

        probes_table.loc[index, mismatch_col] = str(mismatch_probes)
        probes_table.loc[index, probe_length_col] = str(length_probes)

    return probes_table

def calculate_overlaps(probes_table):
    forward_overlap_probe_list = []

    # Forward overlaps.
    for probe in range(0, len(probes_table)):
        overlap_list = []
        for other_probe in range(probe+1, len(probes_table)):
            if (probes_table.iloc[other_probe].name - probes_table.iloc[probe].name) <= 50:
                overlap_list.append(probes_table.iloc[other_probe].name)
        forward_overlap_probe_list.append(overlap_list)

    backward_overlap_probe_list = []

    # Backward overlaps.
    for probe in range(len(probes_table)-1, -1, -1):
        overlap_list = []
        for other_probe in range(probe-1, -1, -1):
            if (probes_table.iloc[probe].name - probes_table.iloc[other_probe].name) <= 50:
                overlap_list.append(probes_table.iloc[other_probe].name)
        backward_overlap_probe_list.insert(0, overlap_list)
        
    probes_table["forward_overlap"] = forward_overlap_probe_list
    probes_table["backward_overlap"] = backward_overlap_probe_list

    return probes_table

def reduce_overlaps(iterate_table, filter_table, direction = "forward", filter = "domain_spanning"):
    if direction == "forward":
        iterate_table = iterate_table
    elif direction == "backward":
        iterate_table = iterate_table[::-1]

    for index, row in iterate_table.iterrows():
        for overlap_probes in row["combined_overlap"]:

            # Check to make sure that we're not checking a overlapping probe that's already been removed.
            if (overlap_probes in filter_table.index) and (index in filter_table.index):
                if filter == "domain_spanning":
                    if filter_table.loc[overlap_probes, "domain_spanning"] < row["domain_spanning"]:
                        print("domain_spanning filter:", index, overlap_probes, row["domain_spanning"], filter_table.loc[overlap_probes, "domain_spanning"])
                        filter_table["combined_overlap"] = [[item for item in inner_list if item != overlap_probes] for inner_list in filter_table["combined_overlap"]]
                        filter_table.drop(index = overlap_probes, inplace = True)

                elif filter == "GC_homopolymer":
                    if filter_table.loc[overlap_probes, "GC_homopolymer"] < row["GC_homopolymer"]:
                        print("GC_homopolymer filter:", index, overlap_probes, row["probe_GC"], filter_table.loc[overlap_probes, "probe_GC"])
                        filter_table["combined_overlap"] = [[item for item in inner_list if item != overlap_probes] for inner_list in filter_table["combined_overlap"]]
                        filter_table.drop(index = overlap_probes, inplace = True)
                
                elif filter == "probe_GC_diff":
                    if row["probe_GC_diff"] == 0 and filter_table.loc[overlap_probes, "probe_GC_diff"] != 0:
                        print("GC Difference filter:", index, overlap_probes, row["probe_GC_diff"], filter_table.loc[overlap_probes, "probe_GC_diff"])
                        filter_table["combined_overlap"] = [[item for item in inner_list if item != overlap_probes] for inner_list in filter_table["combined_overlap"]]
                        filter_table.drop(index = overlap_probes, inplace = True)
                
                elif filter == "GC":
                    if filter_table.loc[overlap_probes, "probe_GC"] < row["probe_GC"]:
                        print("GC filter:", index, overlap_probes, row["probe_GC"], filter_table.loc[overlap_probes, "probe_GC"])
                        filter_table["combined_overlap"] = [[item for item in inner_list if item != overlap_probes] for inner_list in filter_table["combined_overlap"]]
                        filter_table.drop(index = overlap_probes, inplace = True)

                
                # elif final_probes_table.loc[overlap_probes, "probe_homopolymer"] > row["probe_homopolymer"]:
                #     final_probes_table["combined_overlap"] = [[item for item in inner_list if item != overlap_probes] for inner_list in final_probes_table["combined_overlap"]]
                #     print("Repeat filter:", index, overlap_probes, row["probe_homopolymer"], final_probes_table.loc[overlap_probes, "probe_homopolymer"])
                #     final_probes_table.drop(index = overlap_probes, inplace = True)

                # elif final_probes_table.loc[overlap_probes, "probe_GC"] < row["probe_GC"]:
                #     final_probes_table["combined_overlap"] = [[item for item in inner_list if item != overlap_probes] for inner_list in final_probes_table["combined_overlap"]]
                #     print("GC filter:", index, overlap_probes, row["probe_GC"], final_probes_table.loc[overlap_probes, "probe_GC"])
                #     final_probes_table.drop(index = overlap_probes, inplace = True)

                # elif final_probes_table.loc[overlap_probes, "probe_GC"] == row["probe_GC"]:
                #     if final_probes_table.loc[overlap_probes, "probe_homopolymer"] > row["probe_homopolymer"]:
                #         final_probes_table["combined_overlap"] = [[item for item in inner_list if item != overlap_probes] for inner_list in final_probes_table["combined_overlap"]]
                #         print("Repeat filter:", index, overlap_probes, row["probe_homopolymer"], final_probes_table.loc[overlap_probes, "probe_homopolymer"])
                #         final_probes_table.drop(index = overlap_probes, inplace = True)
    
    return filter_table

def dict_compare(d1, d2):
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)

    first_specific = list(d1_keys - shared_keys)
    second_specific = list(d2_keys - shared_keys)
    modified = {o : (d1[o], d2[o]) for o in shared_keys if d1[o] != d2[o]}
    same = set(o for o in shared_keys if d1[o] == d2[o])
    return shared_keys, first_specific, second_specific, modified, same

if __name__ == "__main__":
    # Construct the CAR sequence, where every probe spanning a domain junction is donated with an "N".
    CAR_sense_strand = CAR_sense_strand_creation()

    # The sense strand that we're going to use for the probe creation.
    sense_strand = CAR_sense_strand

    # Create the antisense strand from the sense strand and the reverse sense and antisense strands.
    antisense_strand = strand_creation(sense_strand)
    reverse_sense_strand = sense_strand[::-1]
    reverse_antisense_strand = antisense_strand[::-1]

    # Create the possible probes from the antisense strand.
    probes = probes_from_sequence(reverse_antisense_strand, total_probe_length = 50)
    probes_table = split_probes(probes)
    probes_table["probe_sequence"] = probes_table["probe_sequence"].apply(lambda x: x.replace("N", ""))

    # Calculate some values of the probes.
    probes_table = calculate_probe_values(probes_table)

    # Save the probes table to a csv file.
    probes_table.to_csv("CAR/CAR_all_probes.csv", index=True, header=True)

    # 10X requires a few conditions to be met for the probes to be valid.
    # 1. The 25th nucleotide of the LHS of the probe should be a "t".
    filtered_probes_table = filter_probes(probes_table, filter = "LHS_T_end")

    # 2. The GC content of the LHS and RHS of the probe should be between 44% and 72%.
    filtered_probes_table = filter_probes(filtered_probes_table, filter = "GC_content")

    # 3. Cap the number of homopolymer repeats in the probes. Keep those with less than or equal to 5.
    filtered_probes_table = filter_probes(filtered_probes_table, filter = "homopolymer_repeat")

    # 4. The probes matches to off-target genes should have at least five mismatches in one of the LHS or RHS probes to prevent efficient hybridization.
    # Here we convert the pd probe df to a fasta file, so we can blast it against the human transcriptome.
    probe_fasta(filtered_probes_table)
    
    # 5. Filter the probes that have high off-target homology.
    filtered_probes_table = blast_probes(filtered_probes_table, "CAR/raw_probes/probe/CAR_raw_probe_blast_results.csv")
    print("Blast Filtered Probe Length", len(filtered_probes_table))

    # 6. The probes should not overlap with each other. So, we can remove overlapping probes with worse GC and homopolymer content to get the best possible probes.
    print("Overlap Filters:\n")
    filtered_probes_table = calculate_overlaps(filtered_probes_table)
    
    # Get the union of the forward and backward overlap lists.
    filtered_probes_table["combined_overlap"] = filtered_probes_table["forward_overlap"] + filtered_probes_table["backward_overlap"]
    filtered_probes_table["combined_overlap"] = filtered_probes_table["combined_overlap"].apply(sorted)

    # 7. We now remove the probes that overlap with each other if they have worse GC and homopolymer content.
    # We can do this by iterating through the probes and applying the filters.
    # The backward and forward overlap lists are the same, so we can just iterate through one of them.
    final_filtered_probes_table = copy.deepcopy(filtered_probes_table)
    final_filtered_probes_table = reduce_overlaps(filtered_probes_table, final_filtered_probes_table, direction = "forward" , filter = "domain_spanning")
    final_filtered_probes_table = reduce_overlaps(filtered_probes_table, final_filtered_probes_table, direction = "forward" , filter = "GC")

    # 7. The final probes should now be filtered.
    for index, row in final_filtered_probes_table.iterrows():
        row_LHS_repeats = eval(row["LHS_repeats"])
        row_RHS_repeats = eval(row["RHS_repeats"])

        for overlap_probes in row["combined_overlap"]:
            overlap_LHS_repeat = eval(final_filtered_probes_table.loc[overlap_probes, "LHS_repeats"])
            overlap_RHS_repeat = eval(final_filtered_probes_table.loc[overlap_probes, "RHS_repeats"])

            LHS_shared_keys, LHS_first_specific, LHS_second_specific, LHS_modified, _ = dict_compare(row_LHS_repeats, overlap_LHS_repeat)
            RHS_shared_keys, RHS_first_specific, RHS_second_specific, RHS_modified, _ = dict_compare(row_RHS_repeats, overlap_RHS_repeat)
            print(index, overlap_probes)
            print(LHS_first_specific, LHS_second_specific)
            print(RHS_first_specific, RHS_second_specific)
            print("LHS", LHS_modified)
            print("RHS", RHS_modified)