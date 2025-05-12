#git-branch-test
from Bio import SeqIO
import pandas as pd

def calculate_mutation_rate(reference_seq, sequences):
    sequence_length = len(reference_seq)

    mutation_rates = []
    mutation_details = {i+1: [] for i in range(sequence_length)}  # 初始化字典，键为位置，值为空列表

    for pos in range(sequence_length):
        ref_aa = reference_seq[pos]
        total_count = 0
        mutated_count = 0

        for record in sequences:
            amino_acid = record.seq[pos]
            if amino_acid == '-':
                continue  # 缺失的序列不计入计算
            total_count += 1
            if amino_acid != ref_aa:
                mutated_count += 1
                mutation_details[pos + 1].append((record.id, amino_acid))

        # 计算突变率
        mutation_rate = (mutated_count / total_count) if total_count > 0 else 0
        mutation_rates.append(mutation_rate)

    return mutation_rates, mutation_details

def write_results_to_excel(mutation_rates, mutation_details):
    # 突变率结果表
    df1 = pd.DataFrame({
        'Position': range(1, len(mutation_rates) + 1),
        'Mutation Rate': mutation_rates
    })

    # 突变详情表
    mutation_list = []
    for pos, details in mutation_details.items():
        for seq_id, amino_acid in details:
            mutation_list.append([pos, amino_acid, seq_id])
    df2 = pd.DataFrame(mutation_list, columns=['Position', 'Mutated Amino Acid', 'Sequence Name'])

    with pd.ExcelWriter('mutation_analysis.xlsx') as writer:
        df1.to_excel(writer, sheet_name='Mutation Rates', index=False)
        df2.to_excel(writer, sheet_name='Mutation Details', index=False)

def main(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    reference_seq = sequences[0].seq  # 第一条序列作为参考序列

    mutation_rates, mutation_details = calculate_mutation_rate(reference_seq, sequences[1:])  # 排除第一条参考序列
    write_results_to_excel(mutation_rates, mutation_details)

if __name__ == "__main__":
    fasta_file = "/home/wangyizhen/syp/project/RSVF-statistic/RSVAF_AA_2014-2024.fas"
    main(fasta_file)
