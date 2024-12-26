# This script is copied from home made script 20241220_validate model 07_analysis short long region_03_01.py
# This script is to validate residues of "hard missing", either by "short" or "long" regions or overall.
# This script could also be used to validated residues from "modeled" and "soft missing" group by some modification.





import numpy as np
from collections import Counter
import statistics
from collections import OrderedDict



from collections import Counter

def transform_feature_predicted_list(feature_predicted_list):
    """Transform feature_predicted_list (0 -> -2, 1 -> -1, 2 -> 1)."""
    return [-2 if x == 0 else -1 if x == 1 else 1 for x in feature_predicted_list]


def transform_feature_pLDDT_list(feature_pLDDT_list):
    """Transform feature_pLDDT_list (>= 70 -> 1, < 70 -> -2)."""
    return [1 if x >= 70 else -2 for x in feature_pLDDT_list]


def transform_feature_IUPred_list(feature_IUPred_list):
    """Transform feature_IUPred_list (>= 0.5 -> -2, < 0.5 -> 1)."""
    return [-2 if x >= 0.5 else 1 for x in feature_IUPred_list]

    #adjust threshold for IUPred
    #return [-2 if x >= 0.38 else 1 for x in feature_IUPred_list]


def transform_target_list(feature_list, rule_type="short"):
    """
    Apply the transformation rules on a feature list.
    - For "short": Continuous -2 <= 3 => 'T', else 'N'.
    - For "long": Continuous -2 > 3 => 'T', else 'N'.
    """
    transformed_list = []
    n = len(feature_list)
    for i in range(n):
        if feature_list[i] == 1 or feature_list[i] == -1:
            transformed_list.append("N")
        elif feature_list[i] == -2:
            # Count the window size of continuous -2
            count = 1
            # Look backward
            j = i - 1
            while j >= 0 and feature_list[j] == -2:
                count += 1
                j -= 1
                # Look forward
            j = i + 1
            while j < n and feature_list[j] == -2:
                count += 1
                j += 1
                # Apply the transformation rules


            # region_boundary = 30, then short vs long is at 30 residues
            # if region_boundary is a large number like 3000, then all residues are seen as hard missing
            # (no considering for short or long, every residues of hard missing is counted)
            region_boundary = 30

            if rule_type == "short" and count <= region_boundary:
                transformed_list.append("T")
            elif rule_type == "long" and count > region_boundary:
                transformed_list.append("T")
            else:
                transformed_list.append("N")
    return transformed_list


def compare_transformed_lists(target_list, query_list):
    """
    Compare two transformed lists (target vs query).
    Generate 'TP', 'TN', 'FP', and 'FN' based on the rules.
    """
    comparison = []
    for target, query in zip(target_list, query_list):
        if target == "T" and query == "T":
            comparison.append("TP")
        elif target == "N" and query == "N":
            comparison.append("TN")
        elif target == "N" and query == "T":
            comparison.append("FP")
        elif target == "T" and query == "N":
            comparison.append("FN")
    return comparison

def process_file(input_file, output1, output2, summary_file):
    with open(input_file, "r") as infile, \
         open(output1, "w") as out1, \
         open(output2, "w") as out2, \
         open(summary_file, "w") as summary:

        lines = infile.readlines()
        n = len(lines)

        # Check if number of lines is a multiple of 5
        if n % 6 != 0:
            raise ValueError("Input file must contain lines in multiples of 5.")

        # Initialize counters for summary (to avoid extra memory for accumulation)
        overall_counts = {
            "short_target_vs_short_pLDDT": Counter(),
            "long_target_vs_long_pLDDT": Counter(),
            "short_target_vs_short_IUPred": Counter(),
            "long_target_vs_long_IUPred": Counter(),
            "short_target_vs_short_query": Counter(),
            "long_target_vs_long_query": Counter(),
        }

        short_pLDDT_TP_sequence = []
        short_pLDDT_TN_sequence = []
        short_pLDDT_FP_sequence = []
        short_pLDDT_FN_sequence = []

        short_pLDDT_TP_pLDDT = []
        short_pLDDT_TN_pLDDT = []
        short_pLDDT_FP_pLDDT = []
        short_pLDDT_FN_pLDDT = []

        short_pLDDT_TP_IUPred = []
        short_pLDDT_TN_IUPred = []
        short_pLDDT_FP_IUPred = []
        short_pLDDT_FN_IUPred = []

        long_pLDDT_TP_sequence = []
        long_pLDDT_TN_sequence = []
        long_pLDDT_FP_sequence = []
        long_pLDDT_FN_sequence = []

        long_pLDDT_TP_IUPred = []
        long_pLDDT_TN_IUPred = []
        long_pLDDT_FP_IUPred = []
        long_pLDDT_FN_IUPred = []

        long_pLDDT_TP_pLDDT = []
        long_pLDDT_TN_pLDDT = []
        long_pLDDT_FP_pLDDT = []
        long_pLDDT_FN_pLDDT = []


        short_IUPred_TP_sequence = []
        short_IUPred_TN_sequence = []
        short_IUPred_FP_sequence = []
        short_IUPred_FN_sequence = []

        short_IUPred_TP_pLDDT = []
        short_IUPred_TN_pLDDT = []
        short_IUPred_FP_pLDDT = []
        short_IUPred_FN_pLDDT = []

        short_IUPred_TP_IUPred = []
        short_IUPred_TN_IUPred = []
        short_IUPred_FP_IUPred = []
        short_IUPred_FN_IUPred = []

        long_IUPred_TP_sequence = []
        long_IUPred_TN_sequence = []
        long_IUPred_FP_sequence = []
        long_IUPred_FN_sequence = []

        long_IUPred_TP_pLDDT = []
        long_IUPred_TN_pLDDT = []
        long_IUPred_FP_pLDDT = []
        long_IUPred_FN_pLDDT = []

        long_IUPred_TP_IUPred = []
        long_IUPred_TN_IUPred = []
        long_IUPred_FP_IUPred = []
        long_IUPred_FN_IUPred = []

        short_query_TP_sequence = []
        short_query_TN_sequence = []
        short_query_FP_sequence = []
        short_query_FN_sequence = []

        short_query_TP_pLDDT = []
        short_query_TN_pLDDT = []
        short_query_FP_pLDDT = []
        short_query_FN_pLDDT = []

        short_query_TP_IUPred = []
        short_query_TN_IUPred = []
        short_query_FP_IUPred = []
        short_query_FN_IUPred = []

        long_query_TP_sequence = []
        long_query_TN_sequence = []
        long_query_FP_sequence = []
        long_query_FN_sequence = []

        long_query_TP_pLDDT = []
        long_query_TN_pLDDT = []
        long_query_FP_pLDDT = []
        long_query_FN_pLDDT = []

        long_query_TP_IUPred = []
        long_query_TN_IUPred = []
        long_query_FP_IUPred = []
        long_query_FN_IUPred = []

        # Process groups of 5 lines at a time
        for i in range(0, n, 6):
            try:
                # Read the 5 lines of input
                id_value = lines[i].strip()
                feature_sequence = lines[i+1].strip()
                feature_pLDDT = [float(x) for x in lines[i + 2].strip().split()]
                feature_IUPred = [float(x) for x in lines[i + 3].strip().split()]
                feature_target = [int(x) for x in lines[i + 4].strip().split()]
                feature_predicted = [int(x) for x in lines[i + 5].strip().split()]

                # Transformations
                new_feature_predicted_list = transform_feature_predicted_list(feature_predicted)
                new_feature_pLDDT_list = transform_feature_pLDDT_list(feature_pLDDT)
                new_feature_IUPred_list = transform_feature_IUPred_list(feature_IUPred)

                # Target transformation (short/long)
                trans_short_feature_target_list = transform_target_list(feature_target, rule_type="short")
                trans_long_feature_target_list = transform_target_list(feature_target, rule_type="long")

                # New feature transformation (short/long)
                trans_short_new_feature_predicted_list = transform_target_list(new_feature_predicted_list, rule_type="short")
                trans_long_new_feature_predicted_list = transform_target_list(new_feature_predicted_list, rule_type="long")
                trans_short_new_feature_pLDDT_list = transform_target_list(new_feature_pLDDT_list, rule_type="short")
                trans_long_new_feature_pLDDT_list = transform_target_list(new_feature_pLDDT_list, rule_type="long")
                trans_short_new_feature_IUPred_list = transform_target_list(new_feature_IUPred_list, rule_type="short")
                trans_long_new_feature_IUPred_list = transform_target_list(new_feature_IUPred_list, rule_type="long")

                # Perform comparisons and count results for summary
                def compare_and_count(target_list, query_list, key):
                    comparison = compare_transformed_lists(target_list, query_list)
                    # Update the respective counter in overall_counts
                    overall_counts[key].update(comparison)
                    return comparison

                short_query = compare_and_count(trans_short_feature_target_list, trans_short_new_feature_predicted_list,
                                                "short_target_vs_short_query")
                long_query = compare_and_count(trans_long_feature_target_list, trans_long_new_feature_predicted_list,
                                               "long_target_vs_long_query")
                short_pLDDT = compare_and_count(trans_short_feature_target_list, trans_short_new_feature_pLDDT_list,
                                                "short_target_vs_short_pLDDT")
                long_pLDDT = compare_and_count(trans_long_feature_target_list, trans_long_new_feature_pLDDT_list,
                                               "long_target_vs_long_pLDDT")
                short_IUPred = compare_and_count(trans_short_feature_target_list, trans_short_new_feature_IUPred_list,
                                                 "short_target_vs_short_IUPred")
                long_IUPred = compare_and_count(trans_long_feature_target_list, trans_long_new_feature_IUPred_list,
                                                "long_target_vs_long_IUPred")

                short_pLDDT_TP_sequence.extend(
                    feature_sequence[i] for i in range(len(short_pLDDT)) if short_pLDDT[i] == 'TP')
                short_pLDDT_TN_sequence.extend(feature_sequence[i] for i in range(len(short_pLDDT)) if
                                               short_pLDDT[i] == 'TN')
                short_pLDDT_FP_sequence.extend(feature_sequence[i] for i in range(len(short_pLDDT)) if
                                               short_pLDDT[i] == 'FP')
                short_pLDDT_FN_sequence.extend(feature_sequence[i] for i in range(len(short_pLDDT)) if
                                               short_pLDDT[i] == 'FN')

                short_pLDDT_TP_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_pLDDT)) if
                                            short_pLDDT[i] == 'TP')
                short_pLDDT_TN_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_pLDDT)) if
                                            short_pLDDT[i] == 'TN')
                short_pLDDT_FP_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_pLDDT)) if
                                            short_pLDDT[i] == 'FP')
                short_pLDDT_FN_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_pLDDT)) if
                                            short_pLDDT[i] == 'FN')

                short_pLDDT_TP_IUPred.extend(feature_IUPred[i] for i in range(len(short_pLDDT)) if
                                             short_pLDDT[i] == 'TP')
                short_pLDDT_TN_IUPred.extend(feature_IUPred[i] for i in range(len(short_pLDDT)) if
                                             short_pLDDT[i] == 'TN')
                short_pLDDT_FP_IUPred.extend(feature_IUPred[i] for i in range(len(short_pLDDT)) if
                                             short_pLDDT[i] == 'FP')
                short_pLDDT_FN_IUPred.extend(feature_IUPred[i] for i in range(len(short_pLDDT)) if
                                             short_pLDDT[i] == 'FN')

                long_pLDDT_TP_sequence.extend(
                    feature_sequence[i] for i in range(len(long_pLDDT)) if long_pLDDT[i] == 'TP')
                long_pLDDT_TN_sequence.extend(feature_sequence[i] for i in range(len(long_pLDDT)) if
                                              long_pLDDT[i] == 'TN')
                long_pLDDT_FP_sequence.extend(feature_sequence[i] for i in range(len(long_pLDDT)) if
                                              long_pLDDT[i] == 'FP')
                long_pLDDT_FN_sequence.extend(feature_sequence[i] for i in range(len(long_pLDDT)) if
                                              long_pLDDT[i] == 'FN')

                long_pLDDT_TP_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_pLDDT)) if
                                           long_pLDDT[i] == 'TP')
                long_pLDDT_TN_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_pLDDT)) if
                                           long_pLDDT[i] == 'TN')
                long_pLDDT_FP_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_pLDDT)) if
                                           long_pLDDT[i] == 'FP')
                long_pLDDT_FN_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_pLDDT)) if
                                           long_pLDDT[i] == 'FN')

                long_pLDDT_TP_IUPred.extend(feature_IUPred[i] for i in range(len(long_pLDDT)) if
                                            long_pLDDT[i] == 'TP')
                long_pLDDT_TN_IUPred.extend(feature_IUPred[i] for i in range(len(long_pLDDT)) if
                                            long_pLDDT[i] == 'TN')
                long_pLDDT_FP_IUPred.extend(feature_IUPred[i] for i in range(len(long_pLDDT)) if
                                            long_pLDDT[i] == 'FP')
                long_pLDDT_FN_IUPred.extend(feature_IUPred[i] for i in range(len(long_pLDDT)) if
                                            long_pLDDT[i] == 'FN')
                #############################
                ###################################################################################################

                short_IUPred_TP_sequence.extend(feature_sequence[i] for i in range(len(short_IUPred)) if
                                                short_IUPred[i] == 'TP')
                short_IUPred_TN_sequence.extend(feature_sequence[i] for i in range(len(short_IUPred)) if
                                                short_IUPred[i] == 'TN')
                short_IUPred_FP_sequence.extend(feature_sequence[i] for i in range(len(short_IUPred)) if
                                                short_IUPred[i] == 'FP')
                short_IUPred_FN_sequence.extend(feature_sequence[i] for i in range(len(short_IUPred)) if
                                                short_IUPred[i] == 'FN')

                short_IUPred_TP_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_IUPred)) if
                                             short_IUPred[i] == 'TP')
                short_IUPred_TN_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_IUPred)) if
                                             short_IUPred[i] == 'TN')
                short_IUPred_FP_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_IUPred)) if
                                             short_IUPred[i] == 'FP')
                short_IUPred_FN_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_IUPred)) if
                                             short_IUPred[i] == 'FN')

                short_IUPred_TP_IUPred.extend(feature_IUPred[i] for i in range(len(short_IUPred)) if
                                              short_IUPred[i] == 'TP')
                short_IUPred_TN_IUPred.extend(feature_IUPred[i] for i in range(len(short_IUPred)) if
                                              short_IUPred[i] == 'TN')
                short_IUPred_FP_IUPred.extend(feature_IUPred[i] for i in range(len(short_IUPred)) if
                                              short_IUPred[i] == 'FP')
                short_IUPred_FN_IUPred.extend(feature_IUPred[i] for i in range(len(short_IUPred)) if
                                              short_IUPred[i] == 'FN')

                long_IUPred_TP_sequence.extend(feature_sequence[i] for i in range(len(long_IUPred)) if
                                               long_IUPred[i] == 'TP')
                long_IUPred_TN_sequence.extend(feature_sequence[i] for i in range(len(long_IUPred)) if
                                               long_IUPred[i] == 'TN')
                long_IUPred_FP_sequence.extend(feature_sequence[i] for i in range(len(long_IUPred)) if
                                               long_IUPred[i] == 'FP')
                long_IUPred_FN_sequence.extend(feature_sequence[i] for i in range(len(long_IUPred)) if
                                               long_IUPred[i] == 'FN')

                long_IUPred_TP_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_IUPred)) if
                                            long_IUPred[i] == 'TP')
                long_IUPred_TN_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_IUPred)) if
                                            long_IUPred[i] == 'TN')
                long_IUPred_FP_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_IUPred)) if
                                            long_IUPred[i] == 'FP')
                long_IUPred_FN_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_IUPred)) if
                                            long_IUPred[i] == 'FN')

                long_IUPred_TP_IUPred.extend(feature_IUPred[i] for i in range(len(long_IUPred)) if
                                             long_IUPred[i] == 'TP')
                long_IUPred_TN_IUPred.extend(feature_IUPred[i] for i in range(len(long_IUPred)) if
                                             long_IUPred[i] == 'TN')
                long_IUPred_FP_IUPred.extend(feature_IUPred[i] for i in range(len(long_IUPred)) if
                                             long_IUPred[i] == 'FP')
                long_IUPred_FN_IUPred.extend(feature_IUPred[i] for i in range(len(long_IUPred)) if
                                             long_IUPred[i] == 'FN')

                #####################################################################
                #####################################################################

                short_query_TP_sequence.extend(feature_sequence[i] for i in range(len(short_query)) if
                                               short_query[i] == 'TP')
                short_query_TN_sequence.extend(feature_sequence[i] for i in range(len(short_query)) if
                                               short_query[i] == 'TN')
                short_query_FP_sequence.extend(feature_sequence[i] for i in range(len(short_query)) if
                                               short_query[i] == 'FP')
                short_query_FN_sequence.extend(feature_sequence[i] for i in range(len(short_query)) if
                                               short_query[i] == 'FN')

                short_query_TP_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_query)) if
                                            short_query[i] == 'TP')
                short_query_TN_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_query)) if
                                            short_query[i] == 'TN')
                short_query_FP_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_query)) if
                                            short_query[i] == 'FP')
                short_query_FN_pLDDT.extend(feature_pLDDT[i] for i in range(len(short_query)) if
                                            short_query[i] == 'FN')

                short_query_TP_IUPred.extend(feature_IUPred[i] for i in range(len(short_query)) if
                                             short_query[i] == 'TP')
                short_query_TN_IUPred.extend(feature_IUPred[i] for i in range(len(short_query)) if
                                             short_query[i] == 'TN')
                short_query_FP_IUPred.extend(feature_IUPred[i] for i in range(len(short_query)) if
                                             short_query[i] == 'FP')
                short_query_FN_IUPred.extend(feature_IUPred[i] for i in range(len(short_query)) if
                                             short_query[i] == 'FN')

                long_query_TP_sequence.extend(feature_sequence[i] for i in range(len(long_query)) if
                                              long_query[i] == 'TP')
                long_query_TN_sequence.extend(feature_sequence[i] for i in range(len(long_query)) if
                                              long_query[i] == 'TN')
                long_query_FP_sequence.extend(feature_sequence[i] for i in range(len(long_query)) if
                                              long_query[i] == 'FP')
                long_query_FN_sequence.extend(feature_sequence[i] for i in range(len(long_query)) if
                                              long_query[i] == 'FN')

                long_query_TP_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_query)) if
                                           long_query[i] == 'TP')
                long_query_TN_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_query)) if
                                           long_query[i] == 'TN')
                long_query_FP_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_query)) if
                                           long_query[i] == 'FP')
                long_query_FN_pLDDT.extend(feature_pLDDT[i] for i in range(len(long_query)) if
                                           long_query[i] == 'FN')

                long_query_TP_IUPred.extend(feature_IUPred[i] for i in range(len(long_query)) if
                                            long_query[i] == 'TP')
                long_query_TN_IUPred.extend(feature_IUPred[i] for i in range(len(long_query)) if
                                            long_query[i] == 'TN')
                long_query_FP_IUPred.extend(feature_IUPred[i] for i in range(len(long_query)) if
                                            long_query[i] == 'FP')
                long_query_FN_IUPred.extend(feature_IUPred[i] for i in range(len(long_query)) if
                                            long_query[i] == 'FN')


                ############################################


                #print(short_pLDDT_TP_IUPred)
                # Write to output1.txt
                # out1.write(f"{id_value}\n{lines[i + 1].strip()}\n{lines[i + 2].strip()}\n{lines[i + 3].strip()}\n{lines[i + 4].strip()}\n")
                # out1.write(f"{' '.join(short_pLDDT)}\n{' '.join(long_pLDDT)}\n")
                # out1.write(f"{' '.join(short_IUPred)}\n{' '.join(long_IUPred)}\n")
                # out1.write(f"{' '.join(short_query)}\n{' '.join(long_query)}\n")


                # Write to output2.txt
                # out2.write(f"{id_value}\n{lines[i + 1].strip()}\n{lines[i + 2].strip()}\n{lines[i + 3].strip()}\n{lines[i + 4].strip()}\n")
                # out2.write(f"{' '.join(map(str, trans_short_new_feature_pLDDT_list))}\n{' '.join(map(str, trans_long_new_feature_pLDDT_list))}\n")
                # out2.write(f"{' '.join(map(str, trans_short_new_feature_IUPred_list))}\n{' '.join(map(str, trans_long_new_feature_IUPred_list))}\n")
                # out2.write(f"{' '.join(trans_short_feature_target_list)}\n{' '.join(trans_long_feature_target_list)}\n")
                # out2.write(f"{' '.join(trans_short_new_feature_predicted_list)}\n{' '.join(trans_long_new_feature_predicted_list)}\n")

            except Exception as e:
                print(f"Skipping malformed block at lines {i} to {i+4}: {e}")
                continue

        out1.write(
            f"'short_pLDDT_TP_sequence'\n{''.join(map(str, short_pLDDT_TP_sequence))}\n'long_pLDDT_TP_sequence'\n{''.join(map(str, long_pLDDT_TP_sequence))}\n")
        out1.write(
            f"'short_pLDDT_TN_sequence'\n{''.join(map(str, short_pLDDT_TN_sequence))}\n'long_pLDDT_TN_sequence'\n{''.join(map(str, long_pLDDT_TN_sequence))}\n")
        out1.write(
            f"'short_pLDDT_FP_sequence'\n{''.join(map(str, short_pLDDT_FP_sequence))}\n'long_pLDDT_FP_sequence'\n{''.join(map(str, long_pLDDT_FP_sequence))}\n")
        out1.write(
            f"'short_pLDDT_FN_sequence'\n{''.join(map(str, short_pLDDT_FN_sequence))}\n'long_pLDDT_FN_sequence'\n{''.join(map(str, long_pLDDT_FN_sequence))}\n")

        out1.write(
            f"'short_pLDDT_TP_pLDDT'\n{' '.join(map(str, short_pLDDT_TP_pLDDT))}\n'long_pLDDT_TP_pLDDT'\n{' '.join(map(str, long_pLDDT_TP_pLDDT))}\n")
        out1.write(
            f"'short_pLDDT_TN_pLDDT'\n{' '.join(map(str, short_pLDDT_TN_pLDDT))}\n'long_pLDDT_TN_pLDDT'\n{' '.join(map(str, long_pLDDT_TN_pLDDT))}\n")
        out1.write(
            f"'short_pLDDT_FP_pLDDT'\n{' '.join(map(str, short_pLDDT_FP_pLDDT))}\n'long_pLDDT_FP_pLDDT'\n{' '.join(map(str, long_pLDDT_FP_pLDDT))}\n")
        out1.write(
            f"'short_pLDDT_FN_pLDDT'\n{' '.join(map(str, short_pLDDT_FN_pLDDT))}\n'long_pLDDT_FN_pLDDT'\n{' '.join(map(str, long_pLDDT_FN_pLDDT))}\n")

        out1.write(
            f"'short_pLDDT_TP_IUPred'\n{' '.join(map(str, short_pLDDT_TP_IUPred))}\n'long_pLDDT_TP_IUPred'\n{' '.join(map(str, long_pLDDT_TP_IUPred))}\n")
        out1.write(
            f"'short_pLDDT_TN_IUPred'\n{' '.join(map(str, short_pLDDT_TN_IUPred))}\n'long_pLDDT_TN_IUPred'\n{' '.join(map(str, long_pLDDT_TN_IUPred))}\n")
        out1.write(
            f"'short_pLDDT_FP_IUPred'\n{' '.join(map(str, short_pLDDT_FP_IUPred))}\n'long_pLDDT_FP_IUPred'\n{' '.join(map(str, long_pLDDT_FP_IUPred))}\n")
        out1.write(
            f"'short_pLDDT_FN_IUPred'\n{' '.join(map(str, short_pLDDT_FN_IUPred))}\n'long_pLDDT_FN_IUPred'\n{' '.join(map(str, long_pLDDT_FN_IUPred))}\n")

        out1.write(
            f"'short_IUPred_TP_sequence'\n{''.join(map(str, short_IUPred_TP_sequence))}\n'long_IUPred_TP_sequence'\n{''.join(map(str, long_IUPred_TP_sequence))}\n")
        out1.write(
            f"'short_IUPred_TN_sequence'\n{''.join(map(str, short_IUPred_TN_sequence))}\n'long_IUPred_TN_sequence'\n{''.join(map(str, long_IUPred_TN_sequence))}\n")
        out1.write(
            f"'short_IUPred_FP_sequence'\n{''.join(map(str, short_IUPred_FP_sequence))}\n'long_IUPred_FP_sequence'\n{''.join(map(str, long_IUPred_FP_sequence))}\n")
        out1.write(
            f"'short_IUPred_FN_sequence'\n{''.join(map(str, short_IUPred_FN_sequence))}\n'long_IUPred_FN_sequence'\n{''.join(map(str, long_IUPred_FN_sequence))}\n")

        out1.write(
            f"'short_IUPred_TP_pLDDT'\n{' '.join(map(str, short_IUPred_TP_pLDDT))}\n'long_IUPred_TP_pLDDT'\n{' '.join(map(str, long_IUPred_TP_pLDDT))}\n")
        out1.write(
            f"'short_IUPred_TN_pLDDT'\n{' '.join(map(str, short_IUPred_TN_pLDDT))}\n'long_IUPred_TN_pLDDT'\n{' '.join(map(str, long_IUPred_TN_pLDDT))}\n")
        out1.write(
            f"'short_IUPred_FP_pLDDT'\n{' '.join(map(str, short_IUPred_FP_pLDDT))}\n'long_IUPred_FP_pLDDT'\n{' '.join(map(str, long_IUPred_FP_pLDDT))}\n")
        out1.write(
            f"'short_IUPred_FN_pLDDT'\n{' '.join(map(str, short_IUPred_FN_pLDDT))}\n'long_IUPred_FN_pLDDT'\n{' '.join(map(str, long_IUPred_FN_pLDDT))}\n")

        out1.write(
            f"'short_IUPred_TP_IUPred'\n{' '.join(map(str, short_IUPred_TP_IUPred))}\n'long_IUPred_TP_IUPred'\n{' '.join(map(str, long_IUPred_TP_IUPred))}\n")
        out1.write(
            f"'short_IUPred_TN_IUPred'\n{' '.join(map(str, short_IUPred_TN_IUPred))}\n'long_IUPred_TN_IUPred'\n{' '.join(map(str, long_IUPred_TN_IUPred))}\n")
        out1.write(
            f"'short_IUPred_FP_IUPred'\n{' '.join(map(str, short_IUPred_FP_IUPred))}\n'long_IUPred_FP_IUPred'\n{' '.join(map(str, long_IUPred_FP_IUPred))}\n")
        out1.write(
            f"'short_IUPred_FN_IUPred'\n{' '.join(map(str, short_IUPred_FN_IUPred))}\n'long_IUPred_FN_IUPred'\n{' '.join(map(str, long_IUPred_FN_IUPred))}\n")

        out1.write(
            f"'short_query_TP_sequence'\n{''.join(map(str, short_query_TP_sequence))}\n'long_query_TP_sequence'\n{''.join(map(str, long_query_TP_sequence))}\n")
        out1.write(
            f"'short_query_TN_sequence'\n{''.join(map(str, short_query_TN_sequence))}\n'long_query_TN_sequence'\n{''.join(map(str, long_query_TN_sequence))}\n")
        out1.write(
            f"'short_query_FP_sequence'\n{''.join(map(str, short_query_FP_sequence))}\n'long_query_FP_sequence'\n{''.join(map(str, long_query_FP_sequence))}\n")
        out1.write(
            f"'short_query_FN_sequence'\n{''.join(map(str, short_query_FN_sequence))}\n'long_query_FN_sequence'\n{''.join(map(str, long_query_FN_sequence))}\n")

        out1.write(
            f"'short_query_TP_pLDDT'\n{' '.join(map(str, short_query_TP_pLDDT))}\n'long_query_TP_pLDDT'\n{' '.join(map(str, long_query_TP_pLDDT))}\n")
        out1.write(
            f"'short_query_TN_pLDDT'\n{' '.join(map(str, short_query_TN_pLDDT))}\n'long_query_TN_pLDDT'\n{' '.join(map(str, long_query_TN_pLDDT))}\n")
        out1.write(
            f"'short_query_FP_pLDDT'\n{' '.join(map(str, short_query_FP_pLDDT))}\n'long_query_FP_pLDDT'\n{' '.join(map(str, long_query_FP_pLDDT))}\n")
        out1.write(
            f"'short_query_FN_pLDDT'\n{' '.join(map(str, short_query_FN_pLDDT))}\n'long_query_FN_pLDDT'\n{' '.join(map(str, long_query_FN_pLDDT))}\n")

        out1.write(
            f"'short_query_TP_IUPred'\n{' '.join(map(str, short_query_TP_IUPred))}\n'long_query_TP_IUPred'\n{' '.join(map(str, long_query_TP_IUPred))}\n")
        out1.write(
            f"'short_query_TN_IUPred'\n{' '.join(map(str, short_query_TN_IUPred))}\n'long_query_TN_IUPred'\n{' '.join(map(str, long_query_TN_IUPred))}\n")
        out1.write(
            f"'short_query_FP_IUPred'\n{' '.join(map(str, short_query_FP_IUPred))}\n'long_query_FP_IUPred'\n{' '.join(map(str, long_query_FP_IUPred))}\n")
        out1.write(
            f"'short_query_FN_IUPred'\n{' '.join(map(str, short_query_FN_IUPred))}\n'long_query_FN_IUPred'\n{' '.join(map(str, long_query_FN_IUPred))}\n")

        list_names = [
            "short_pLDDT_TP_sequence",
            "short_pLDDT_TN_sequence",
            "short_pLDDT_FP_sequence",
            "short_pLDDT_FN_sequence",

            "short_pLDDT_TP_pLDDT",
            "short_pLDDT_TN_pLDDT",
            "short_pLDDT_FP_pLDDT",
            "short_pLDDT_FN_pLDDT",

            "short_pLDDT_TP_IUPred",
            "short_pLDDT_TN_IUPred",
            "short_pLDDT_FP_IUPred",
            "short_pLDDT_FN_IUPred",

            "long_pLDDT_TP_sequence",
            "long_pLDDT_TN_sequence",
            "long_pLDDT_FP_sequence",
            "long_pLDDT_FN_sequence",

            "long_pLDDT_TP_IUPred",
            "long_pLDDT_TN_IUPred",
            "long_pLDDT_FP_IUPred",
            "long_pLDDT_FN_IUPred",

            "long_pLDDT_TP_pLDDT",
            "long_pLDDT_TN_pLDDT",
            "long_pLDDT_FP_pLDDT",
            "long_pLDDT_FN_pLDDT",

            "short_IUPred_TP_sequence",
            "short_IUPred_TN_sequence",
            "short_IUPred_FP_sequence",
            "short_IUPred_FN_sequence",

            "short_IUPred_TP_pLDDT",
            "short_IUPred_TN_pLDDT",
            "short_IUPred_FP_pLDDT",
            "short_IUPred_FN_pLDDT",

            "short_IUPred_TP_IUPred",
            "short_IUPred_TN_IUPred",
            "short_IUPred_FP_IUPred",
            "short_IUPred_FN_IUPred",

            "long_IUPred_TP_sequence",
            "long_IUPred_TN_sequence",
            "long_IUPred_FP_sequence",
            "long_IUPred_FN_sequence",

            "long_IUPred_TP_pLDDT",
            "long_IUPred_TN_pLDDT",
            "long_IUPred_FP_pLDDT",
            "long_IUPred_FN_pLDDT",

            "long_IUPred_TP_IUPred",
            "long_IUPred_TN_IUPred",
            "long_IUPred_FP_IUPred",
            "long_IUPred_FN_IUPred",

            "short_query_TP_sequence",
            "short_query_TN_sequence",
            "short_query_FP_sequence",
            "short_query_FN_sequence",

            "short_query_TP_pLDDT",
            "short_query_TN_pLDDT",
            "short_query_FP_pLDDT",
            "short_query_FN_pLDDT",

            "short_query_TP_IUPred",
            "short_query_TN_IUPred",
            "short_query_FP_IUPred",
            "short_query_FN_IUPred",

            "long_query_TP_sequence",
            "long_query_TN_sequence",
            "long_query_FP_sequence",
            "long_query_FN_sequence",

            "long_query_TP_pLDDT",
            "long_query_TN_pLDDT",
            "long_query_FP_pLDDT",
            "long_query_FN_pLDDT",

            "long_query_TP_IUPred",
            "long_query_TN_IUPred",
            "long_query_FP_IUPred",
            "long_query_FN_IUPred"
        ]




        tracking_characters = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F',
                               'Y', 'W']

        # Open the output file
        for name in list_names:
            mean_name = f"{name}_mean"
            std_name = f"{name}_std"

            # Check if the list ends with "_sequence"
            if name.endswith("_sequence"):
                # Count occurrences of each character we are tracking
                counts = {char: 0 for char in tracking_characters}
                for char in eval(name):
                    if char in counts:
                        counts[char] += 1

                        # Combine counts into the specified output format
                counts_str = '\t'.join(str(counts[char]) for char in tracking_characters)
                character_str = '\t'.join(tracking_characters)

                # Write out the results in the desired format
                out2.write(f"{name}: {character_str}\n{counts_str}\n")
            else:
                # Calculate mean and standard deviation for non-sequence lists
                mean_value = np.mean(eval(name))
                std_value = np.std(eval(name))
                out2.write(f"'{mean_name}'\t'{std_name}'\n")
                out2.write(f"{mean_value}\t{std_value}\n")




                #
        # print("Mean:", short_query_TP_IUPred_mean)
        # print("Standard Deviation:", short_query_TP_IUPred_std)


        # Write summary1.txt
        for key, count in overall_counts.items():
            summary.write(f"{key} TP is :\t{count['TP']}\n")
            summary.write(f"{key} TN is :\t{count['TN']}\n")
            summary.write(f"{key} FP is :\t{count['FP']}\n")
            summary.write(f"{key} FN is :\t{count['FN']}\n\n")



            # Example usage



basic_path_example = 'D:\\Example\\04_validate_results\\'
#base_path = 'C:\\Users\\42572\\Desktop\\data_figures\\training\\best processed\\'
base_path = basic_path_example

base_input = 'Validation_result_Xray before 2022_output2_balanced_model_256_1500_e2_Xray behind 2022_output2_AAplus'
base_output1 = 'outputs\\' + base_input + '_output1_short_long'
base_output2 = 'outputs\\' + base_input + '_output2_short_long'
base_summar1 = 'outputs\\' + base_input + '_summary_short_long'

input = base_path + 'inputs\\' + base_input + '.txt'
output1 = base_path  + base_output1 + '.txt'
output2 = base_path  + base_output2 + '.txt'
summary1 = base_path + base_summar1 + '.txt'



process_file(input, output1, output2, summary1)