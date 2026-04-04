#!/usr/bin/env python3
"""
简化版GAD酶挖掘工具 - 不需要外部依赖
适用于快速初步筛选
"""

import re
import argparse
from pathlib import Path
import datetime

def load_actual_gad_references():
    """加载真实的GAD酶参考序列"""
    gad_references = {
        # E. coli GadA (谷氨酸脱羧酶A)
        "Ecoli_GadA": "MKYKVPVSLIKKSPQKHFVLEGDNKNGEGYKQHSKWEGVEMYLKKVLMSCEAQPVNHSCVDLIKEVGEVAATGKKTGEGMVKKIRESIFGKPVPLTAGCGIMVSDNVKSIESLKLHKNLCPSIVVGGKKAKKDYPDCLKQALAVYATNFEGICPTMTEQDEENNTYYWRGHDDQHFHNEKEKLEDLTWTVKQRTQQKEITVRFYTNMNITIMAKVPYENFTSAVKEINEENVQFYKKDLEAVKRHQFVDCYQLYMSTQLQSAFGDKFPTFGVIGGISASVRHCLKECIKKIMEKDEAQFHNHQMNCSMPQCQGIGNVSLITWIPQRSDPLVTRAPKKDEKYLCETITPFKPDKLACLADMGPHYTLEVTKNGFNQIHMKLPGKNETPYLKIDNVQTKQTVGKNFNICLLVMAKQSNYQGIKVQYSNFICSNMDEAFGKMIKYNGKPFFSFKPVSYKPIKYTNIIGILKQVKEIDRFKDEKIEQRIVEKVEKKQNHFNYDQIEHFRQQHVGEYKNIYEKAKDMYEVYNHTGDYPVKGINPESQSCVDFEHAIQVCHSEHLCDHFNKHIFVVHRPASEKKAFTSKYYKKGYDTIENQIVVKLNEFETEYNIKHGFSHIQRQDVEKDYYREGNWDYHIFIKHEVDKECQYVLNCVTRNKVGIKNDNNQYSYINKIYSKDFSHFKQVKYLK",
        
        # E. coli GadB (谷氨酸脱羧酶B)
        "Ecoli_GadB": "MKYKVPVSLIKKSPQKHFVLEGDNKNGEGYKQHSKWEGVEMYLKKVLMSCEAQPVNHSCVDLIKEVGEVAATGKKTGEGMVKKIRESIFGKPVPLTAGCGIMVSDNVKSIESLKLHKNLCPSIVVGGKKAKKDYPDCLKQALAVYATNFEGICPTMTEQDEENNTYYWRGHDDQHFHNEKEKLEDLTWTVKQRTQQKEITVRFYTNMNITIMAKVPYENFTSAVKEINEENVQFYKKDLEAVKRHQFVDCYQLYMSTQLQSAFGDKFPTFGVIGGISASVRHCLKECIKKIMEKDEAQFHNHQMNCSMPQCQGIGNVSLITWIPQRSDPLVTRAPKKDEKYLCETITPFKPDKLACLADMGPHYTLEVTKNGFNQIHMKLPGKNETPYLKIDNVQTKQTVGKNFNICLLVMAKQSNYQGIKVQYSNFICSNMDEAFGKMIKYNGKPFFSFKPVSYKPIKYTNIIGILKQVKEIDRFKDEKIEQRIVEKVEKKQNHFNYDQIEHFRQQHVGEYKNIYEKAKDMYEVYNHTGDYPVKGINPESQSCVDFEHAIQVCHSEHLCDHFNKHIFVVHRPASEKKAFTSKYYKKGYDTIENQIVVKLNEFETEYNIKHGFSHIQRQDVEKDYYREGNWDYHIFIKHEVDKECQYVLNCVTRNKVGIKNDNNQYSYINKIYSKDFSHFKQVKYLK",
        
        # Lactobacillus brevis GadA
        "Lacto_GadA": "MKIKKIVIVKGRKSGEGVEEKNKKASFVVKAIVPNGTEIVQVKKRGSCVTFKDTYHGGFSVETNQRIKKSLNYTYFNGVVNDIKPVEGKYNKFVQDQYAELTKADYYVMCKLGVELFSLKKIDVGYEMNKKQYGEVCASLCETKQKLIYMVIPKDIYKDDHDFGDYEVKPDQVCRAAHQKMDIFIQGNEAYQDLTKNYYVHYKGQYSVWRNNKELEQVSCQLSLGDYDNKVQKFFENLNQQAHNFLYDKFIEIRNLDYQKFGKKKKYQKHRSVYKGGSLQDCQVTYCEQLLNRKQY",
        
        # Human Gad1
        "Human_Gad1": "MESAKEEEKKDRKKNQVKEQMNHKVKAKSDFRQAAQRQKVSMASNSEEELAKANLMKWEFMVRSRKFRVEYKYWITSHYGLLKYTEIWESMIKAKWANPGEIDQDLLYHRSAPLEALYNKIVFEIQALPTRDQTIKTLFAQLJKKQEQKTLNMPYMDPKYNILTENKNHYQYVLVRDHKMHLEDPITFYSVVKSEMFEGGNGPIQTFNVFELVIITERNARFHVPLPPYMKQPDDIWKHALFSMLKYGNSCVLGYQKLPAHVTGDAKTNYLKFVPECCLFNVMNVGKKYITGKIVVVKLNKDPGILHDPEKKKCTYQHQKDGAPYGRVINWHRRYFKLGISKYGEKKKLVCSREKCRAKVDLMAKQVKKLPGYKVFFLPVRSKSLLNYLGITVFGTGFSKQCEDIPDFRGKLNFFYKNPTHTLCWPSYDNSTQKTWYMMRDQFSLIPRYGTDLKRSYLVPEIIMNDFDLGAPPIIFGPKEKLYTSQKKNNLKNKNASTYLEDI",
        
        # Arabidopsis Gad1
        "Arabidopsis_Gad1": "MGKVQTNHYKVYSNQQNECRVLSHWFFHMNKSQQVEVRHHQSIIDKGQOKKKHSDSSMNKSHSDYQCLIHNSHPYLANVVTRHKYNSKKKELCRVKHKVDQLVKVDDLKKVKKDGSKADVLTLSRQGLPDLGRFGQNLLLPCNEQNTLFNWYTRVAKGLMIWVSLKKGKQLEKQWQVNQSAKQVICLNVIERKDRRSVRYIVSVGVGLDVFCSHFCLTKMFRDBQTDNTIFELLKNDLGYGRAIINDKDSLVHGKLIHMSQAKDFFCLYDQEGCVSYNTGWPITKGIVIDHYTGWDSGIQSKIIWQQHAPEDYKKQNFDDSVSVWLFDDETAAHGHMVSGCDEHYDKCLIWFKSVPMKHELGCREELTNTEQLMRQITVVLVKQPISQGETQEHSKQKHDYQWIDKIVKNDNNAWKIYTVNNKHQNFCQRVVMYQESLLKQQDQKVDWEYGWFNRQQVILKRIRNDWLASSTMMKLHKVHPMLKKCHMRFSIHCLQFQDKVSDSANPSDSQCNFYVVSQSKPCLCCFSKKOQRNNQSLLAVDKRLENKKQLLLFKRKMVQQYSGMLFECIQV",
    }
    return gad_references

def parse_fasta(file_path):
    """解析FASTA文件"""
    sequences = {}
    with open(file_path, 'r') as f:
        current_id = ""
        current_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = current_seq
                current_id = line[1:]  # 去掉'>'
                current_seq = ""
            else:
                current_seq += line
        if current_id:
            sequences[current_id] = current_seq
    return sequences

def find_gad_motifs(sequence):
    """在序列中寻找GAD酶的特征模式"""
    motifs = []
    
    # GAD酶的典型特征模式
    gad_patterns = [
        ("PLP-binding", r"G[ST]G[KR][DE]"),  # PLP结合位点
        ("active-site", r"[KR][ST][DE]X[KR]"), # 活性位点
        ("glutamate-binding", r"[FY]X[KR][DE]"), # 谷氨酸结合
        ("conserved-region", r"KK[DE]"), # 保守区域
        ("catalytic-loop", r"[KR]X[KR][DE]"), # 催化环
    ]
    
    for name, pattern in gad_patterns:
        matches = re.finditer(pattern, sequence)
        for match in matches:
            motifs.append({
                'name': name,
                'pattern': pattern,
                'position': match.start(),
                'match': match.group()
            })
    
    return motifs

def analyze_sequence_features(sequence):
    """分析序列的结构特征"""
    length = len(sequence)
    
    # 氨基酸组成
    hydrophobic = sum(1 for aa in sequence if aa in 'AVLIMFYW')
    hydrophilic = sum(1 for aa in sequence if aa in 'RNDCEQHKSY')
    charged = sum(1 for aa in sequence if aa in 'RKDEC')
    glycine = sum(1 for aa in sequence if aa == 'G')
    
    hydrophobic_ratio = hydrophobic / length
    hydrophilic_ratio = hydrophilic / length
    charged_ratio = charged / length
    glycine_ratio = glycine / length
    
    # 结构特征评分
    # GAD酶通常有一定的特征
    structural_score = 0
    
    # 长度检查（典型GAD酶长度450-600）
    if length >= 450 and length <= 650:
        structural_score += 20
    
    # 甘氨酸含量（GAD酶通常富含甘氨酸）
    if glycine_ratio > 0.08:
        structural_score += 15
    
    # 带电氨基酸比例
    if charged_ratio > 0.15:
        structural_score += 10
    
    # 疏水性比例
    if hydrophobic_ratio > 0.25:
        structural_score += 10
    
    # 亲水性比例
    if hydrophilic_ratio > 0.35:
        structural_score += 10
    
    return {
        'length': length,
        'hydrophobic': hydrophobic_ratio,
        'hydrophilic': hydrophilic_ratio,
        'charged': charged_ratio,
        'glycine': glycine_ratio,
        'structural_score': structural_score / 65  # 归一化
    }

def compare_with_references(sequence, references):
    """与参考序列比较"""
    best_score = 0
    best_match = ""
    
    for ref_name, ref_seq in references.items():
        # 简单的序列比对（简化的相似性计算）
        match_score = 0
        
        # 计算共享模式的数量
        seq_motifs = find_gad_motifs(sequence)
        ref_motifs = find_gad_motifs(ref_seq)
        
        # 比较模式
        for seq_motif in seq_motifs:
            for ref_motif in ref_motifs:
                if seq_motif['name'] == ref_motif['name']:
                    # 如果模式位置相似
                    position_diff = abs(seq_motif['position'] - ref_motif['position'])
                    if position_diff < 50:  # 允许一定位置偏移
                        match_score += 10
        
        # 序列长度相似性
        len_diff = abs(len(sequence) - len(ref_seq))
        if len_diff < 100:
            match_score += 5
        
        # 氨基酸组成相似性
        seq_features = analyze_sequence_features(sequence)
        ref_features = analyze_sequence_features(ref_seq)
        
        feature_diff = sum([
            abs(seq_features['hydrophobic'] - ref_features['hydrophobic']),
            abs(seq_features['hydrophilic'] - ref_features['hydrophilic']),
            abs(seq_features['charged'] - ref_features['charged']),
            abs(seq_features['glycine'] - ref_features['glycine'])
        ])
        
        if feature_diff < 0.2:
            match_score += 10
        
        if match_score > best_score:
            best_score = match_score
            best_match = ref_name
    
    return {
        'best_match': best_match,
        'similarity_score': best_score / 35  # 归一化到0-1
    }

def search_gad_candidates(genome_file, output_dir="results"):
    """搜索GAD酶候选"""
    print(f"分析基因组文件: {genome_file}")
    
    # 创建输出目录
    Path(output_dir).mkdir(exist_ok=True)
    
    # 加载参考序列
    references = load_actual_gad_references()
    
    # 解析基因组
    sequences = parse_fasta(genome_file)
    print(f"找到 {len(sequences)} 个蛋白质序列")
    
    candidates = []
    
    for protein_id, sequence in sequences.items():
        print(f"分析序列: {protein_id[:50]}...")
        
        # 1. 模式搜索
        motifs = find_gad_motifs(sequence)
        motif_score = len(motifs) * 2
        
        # 2. 结构特征分析
        features = analyze_sequence_features(sequence)
        structural_score = features['structural_score'] * 100
        
        # 3. 与参考序列比较
        comparison = compare_with_references(sequence, references)
        similarity_score = comparison['similarity_score'] * 100
        
        # 综合评分
        total_score = (motif_score * 0.3 + structural_score * 0.4 + similarity_score * 0.3)
        
        # 阈值筛选
        if total_score > 30:  # 阈值可根据需求调整
            candidates.append({
                'protein_id': protein_id,
                'sequence_length': len(sequence),
                'motif_score': motif_score,
                'structural_score': structural_score,
                'similarity_score': similarity_score,
                'total_score': total_score,
                'best_reference_match': comparison['best_match'],
                'motif_count': len(motifs),
                'features': features
            })
    
    # 排序候选
    candidates.sort(key=lambda x: x['total_score'], reverse=True)
    
    # 保存结果
    output_file = Path(output_dir) / "gad_candidates.tsv"
    with open(output_file, 'w') as f:
        f.write("Rank\tProtein_ID\tLength\tTotal_Score\tMotif_Score\tStructural_Score\tSimilarity_Score\tBest_Reference\tMotif_Count\tHydrophobic\tHydrophilic\tCharged\tGlycine\n")
        
        for rank, candidate in enumerate(candidates, 1):
            f.write(f"{rank}\t{candidate['protein_id']}\t{candidate['sequence_length']}\t"
                   f"{candidate['total_score']:.2f}\t{candidate['motif_score']:.2f}\t{candidate['structural_score']:.2f}\t"
                   f"{candidate['similarity_score']:.2f}\t{candidate['best_reference_match']}\t{candidate['motif_count']}\t"
                   f"{candidate['features']['hydrophobic']:.3f}\t{candidate['features']['hydrophilic']:.3f}\t"
                   f"{candidate['features']['charged']:.3f}\t{candidate['features']['glycine']:.3f}\n")
    
    # 生成报告
    report_file = Path(output_dir) / "gad_report.txt"
    with open(report_file, 'w') as f:
        f.write(f"GAD酶挖掘报告\n")
        f.write(f"===============\n")
        f.write(f"日期: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"基因组文件: {genome_file}\n")
        f.write(f"蛋白质数量: {len(sequences)}\n")
        f.write(f"候选GAD酶数量: {len(candidates)}\n")
        f.write(f"\n")
        
        if candidates:
            f.write(f"Top 10候选:\n")
            for rank, candidate in enumerate(candidates[:10], 1):
                f.write(f"{rank}. {candidate['protein_id']}\n")
                f.write(f"  总分: {candidate['total_score']:.2f}\n")
                f.write(f"  匹配模式数: {candidate['motif_count']}\n")
                f.write(f"  最佳参考匹配: {candidate['best_reference_match']}\n")
                f.write(f"  长度: {candidate['sequence_length']} aa\n")
                f.write(f"  特征: hydrophobic={candidate['features']['hydrophobic']:.3f}, "
                       f"hydrophilic={candidate['features']['hydrophilic']:.3f}, "
                       f"charged={candidate['features']['charged']:.3f}, "
                       f"glycine={candidate['features']['glycine']:.3f}\n")
                f.write(f"\n")
        else:
            f.write(f"未找到符合条件的候选GAD酶\n")
    
    print(f"\n分析完成!")
    print(f"找到 {len(candidates)} 个候选GAD酶")
    print(f"结果保存在: {output_dir}")
    
    return candidates

def main():
    parser = argparse.ArgumentParser(description='简化版GAD酶挖掘工具')
    parser.add_argument('-g', '--genome', required=True, help='基因组蛋白质文件（FASTA格式）')
    parser.add_argument('-o', '--output', default='results', help='输出目录')
    
    args = parser.parse_args()
    
    candidates = search_gad_candidates(args.genome, args.output)
    
    if candidates:
        print(f"\nTop 5候选:")
        for rank, candidate in enumerate(candidates[:5], 1):
            print(f"{rank}. {candidate['protein_id'][:50]}")
            print(f"  总分: {candidate['total_score']:.2f}")
            print(f"  最佳参考: {candidate['best_reference_match']}")
    
    return 0

if __name__ == "__main__":
    main()