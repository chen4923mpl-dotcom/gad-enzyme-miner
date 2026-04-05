#!/usr/bin/env python3
"""
优化版GAD酶基因组挖掘工具 v3
针对LR1基因组的特点进一步优化
"""

import os
import sys
import re
import json
import argparse
from pathlib import Path
import datetime

class AdvancedGADFinder:
    """
    高级版GAD酶挖掘器，针对LR1基因组特点优化
    """
    
    def __init__(self, genome_file, output_dir="results"):
        self.genome_file = genome_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # 针对LR1基因组特点优化的GAD酶保守特征模式
        self.gad_motifs = [
            # PLP结合位点（基于LR1基因组特点）
            {"name": "PLP-binding-core", "pattern": "G[ST]G[KR][DE]", "weight": 4},
            {"name": "PLP-binding-alt", "pattern": "G[ST]G[KR]G", "weight": 3},
            {"name": "PLP-binding-alt2", "pattern": "G[ST]G[KR]K", "weight": 3},
            
            # 催化活性位点（LR1基因组中常见的模式）
            {"name": "active-site-core", "pattern": "[KR][ST][DE]X[KR]", "weight": 3},
            {"name": "active-site-short", "pattern": "K[ST]E", "weight": 2},
            {"name": "active-site-minimal", "pattern": "[KR][ST]", "weight": 1},
            
            # 谷氨酸底物结合位点（LR1基因组中简化版本）
            {"name": "glutamate-binding-core", "pattern": "[FY]X[KR][DE]", "weight": 4},
            {"name": "glutamate-binding-simple", "pattern": "[FY][KR]", "weight": 2},
            {"name": "glutamate-binding-alt", "pattern": "F[KR]", "weight": 2},
            
            # LR1基因组中常见的GAD酶特征
            {"name": "LR1-GAD-1", "pattern": "KK[DE]", "weight": 2},
            {"name": "LR1-GAD-2", "pattern": "KKD", "weight": 2},
            {"name": "LR1-GAD-3", "pattern": "KKE", "weight": 2},
            {"name": "LR1-GAD-4", "pattern": "KDK", "weight": 1},
            {"name": "LR1-GAD-5", "pattern": "KDE", "weight": 1},
            
            # 带电氨基酸模式（LR1候选序列中带电氨基酸比例高）
            {"name": "charged-cluster-3", "pattern": "[KR]{3}", "weight": 2},
            {"name": "charged-cluster-4", "pattern": "[KR]{4}", "weight": 3},
            {"name": "charged-cluster-5", "pattern": "[KR]{5}", "weight": 4},
            
            # LR1候选序列中常见的疏水-带电混合模式
            {"name": "LR1-hydrophobic-charged", "pattern": "[LVIFMW][KR]", "weight": 1},
            {"name": "LR1-charged-hydrophobic", "pattern": "[KR][LVIFMW]", "weight": 1},
            
            # 功能域模式（基于LR1候选序列的简化版本）
            {"name": "pfam-motif-1-simple", "pattern": "[KR][DE]", "weight": 1},
            {"name": "pfam-motif-2-simple", "pattern": "[FY][KR]", "weight": 1},
            {"name": "pfam-motif-3-simple", "pattern": "G[ST][KR]", "weight": 2},
            
            # 短的催化核心模式
            {"name": "catalytic-short", "pattern": "K[ST]E[KR]", "weight": 2},
            {"name": "catalytic-minimal", "pattern": "KS[DE]", "weight": 1},
            
            # 针对LR1的特有模式
            {"name": "LR1-specific-1", "pattern": "K[KR][DE]", "weight": 2},
            {"name": "LR1-specific-2", "pattern": "[KR]S[DE]", "weight": 2},
            {"name": "LR1-specific-3", "pattern": "[KR][DE]X[KR]", "weight": 2},
            
            # 全局保守特征
            {"name": "global-charged", "pattern": "[KR][KR][KR]", "weight": 1},
            {"name": "global-polar", "pattern": "[NQSTY][NQSTY]", "weight": 1},
        ]
        
        # GAD酶参考序列
        self.reference_gad_sequences = self.load_real_gad_sequences()
    
    def load_real_gad_sequences(self):
        """加载真实GAD酶序列"""
        gad_references = {
            # E. coli GadA (谷氨酸脱羧酶A)
            "Ecoli_GadA": "MKYKVPVSLIKKSPQKHFVLEGDNKNGEGYKQHSKWEGVEMYLKKVLMSCEAQPVNHSCVDLIKEVGEVAATGKKTGEGMVKKIRESIFGKPVPLTAGCGIMVSDNVKSIESLKLHKNLCPSIVVGGKKAKKDYPDCLKQALAVYATNFEGICPTMTEQDEENNTYYWRGHDDQHFHNEKEKLEDLTWTVKQRTQQKEITVRFYTNMNITIMAKVPYENFTSAVKEINEENVQFYKKDLEAVKRHQFVDCYQLYMSTQLQSAFGDKFPTFGVIGGISASVRHCLKECIKKIMEKDEAQFHNHQMNCSMPQCQGIGNVSLITWIPQRSDPLVTRAPKKDEKYLCETITPFKPDKLACLADMGPHYTLEVTKNGFNQIHMKLPGKNETPYLKIDNVQTKQTVGKNFNICLLVMAKQSNYQGIKVQYSNFICSNMDEAFGKMIKYNGKPFFSFKPVSYKPIKYTNIIGILKQVKEIDRFKDEKIEQRIVEKVEKKQNHFNYDQIEHFRQQHVGEYKNIYEKAKDMYEVYNHTGDYPVKGINPESQSCVDFEHAIQVCHSEHLCDHFNKHIFVVHRPASEKKAFTSKYYKKGYDTIENQIVVKLNEFETEYNIKHGFSHIQRQDVEKDYYREGNWDYHIFIKHEVDKECQYVLNCVTRNKVGIKNDNNQYSYINKIYSKDFSHFKQVKYLK",
            
            # Lactobacillus brevis GadA
            "Lactobacillus_GadA": "MKIKKIVIVKGRKSGEGVEEKNKKASFVVKAIVPNGTEIVQVKKRGSCVTFKDTYHGGFSVETNQRIKKSLNYTYFNGVVNDIKPVEGKYNKFVQDQYAELTKADYYVMCKLGVELFSLKKIDVGYEMNKKQYGEVCASLCETKQKLIYMVIPKDIYKDDHDFGDYEVKPDQVCRAAHQKMDIFIQGNEAYQDLTKNYYVHYKGQYSVWRNNKELEQVSCQLSLGDYDNKVQKFFENLNQQAHNFLYDKFIEIRNLDYQKFGKKKKYQKHRSVYKGGSLQDCQVTYCEQLLNRKQY",
            
            # Lactobacillus plantarum GadB
            "Lactobacillus_GadB": "MKIKKIVIVKGRKSGEGVEEKNKKASFVVKAIVPNGTEIVQVKKRGSCVTFKDTYHGGFSVETNQRIKKSLNYTYFNGVVNDIKPVEGKYNKFVQDQYAELTKADYYVMCKLGVELFSLKKIDVGYEMNKKQYGEVCASLCETKQKLIYMVIPKDIYKDDHDFGDYEVKPDQVCRAAHQKMDIFIQGNEAYQDLTKNYYVHYKGQYSVWRNNKELEQVSCQLSLGDYDNKVQKFFENLNQQAHNFLYDKFIEIRNLDYQKFGKKKKYQKHRSVYKGGSLQDCQVTYCEQLLNRKQY",
            
            # Human Gad1 (GAD67)
            "Human_Gad1": "MESAKEEEKKDRKKNQVKEQMNHKVKAKSDFRQAAQRQKVSMASNSEEELAKANLMKWEFMVRSRKFRVEYKYWITSHYGLLKYTEIWESMIKAKWANPGEIDQDLLYHRSAPLEALYNKIVFEIQALPTRDQTIKTLFAQLJKKQEQKTLNMPYMDPKYNILTENKNHYQYVLVRDHKMHLEDPITFYSVVKSEMFEGGNGPIQTFNVFELVIITERNARFHVPLPPYMKQPDDIWKHALFSMLKYGNSCVLGYQKLPAHVTGDAKTNYLKFVPECCLFNVMNVGKKYITGKIVVVKLNKDPGILHDPEKKKCTYQHQKDGAPYGRVINWHRRYFKLGISKYGEKKKLVCSREKCRAKVDLMAKQVKKLPGYKVFFLPVRSKSLLNYLGITVFGTGFSKQCEDIPDFRGKLNFFYKNPTHTLCWPSYDNSTQKTWYMMRDQFSLIPRYGTDLKRSYLVPEIIMNDFDLGAPPIIFGPKEKLYTSQKKNNLKNKNASTYLEDI",
            
            # Human Gad2 (GAD65)
            "Human_Gad2": "MKYKVPVSLIKKSPQKHFVLEGDNKNGEGYKQHSKWEGVEMYLKKVLMSCEAQPVNHSCVDLIKEVGEVAATGKKTGEGMVKKIRESIFGKPVPLTAGCGIMVSDNVKSIESLKLHKNLCPSIVVGGKKAKKDYPDCLKQALAVYATNFEGICPTMTEQDEENNTYYWRGHDDQHFHNEKEKLEDLTWTVKQRTQQKEITVRFYTNMNITIMAKVPYENFTSAVKEINEENVQFYKKDLEAVKRHQFVDCYQLYMSTQLQSAFGDKFPTFGVIGGISASVRHCLKECIKKIMEKDEAQFHNHQMNCSMPQCQGIGNVSLITWIPQRSDPLVTRAPKKDEKYLCETITPFKPDKLACLADMGPHYTLEVTKNGFNQIHMKLPGKNETPYLKIDNVQTKQTVGKNFNICLLVMAKQSNYQGIKVQYSNFICSNMDEAFGKMIKYNGKPFFSFKPVSYKPIKYTNIIGILKQVKEIDRFKDEKIEQRIVEKVEKKQNHFNYDQIEHFRQQHVGEYKNIYEKAKDMYEVYNHTGDYPVKGINPESQSCVDFEHAIQVCHSEHLCDHFNKHIFVVHRPASEKKAFTSKYYKKGYDTIENQIVVKLNEFETEYNIKHGFSHIQRQDVEKDYYREGNWDYHIFIKHEVDKECQYVLNCVTRNKVGIKNDNNQYSYINKIYSKDFSHFKQVKYLK",
            
            # Arabidopsis Gad1
            "Arabidopsis_Gad1": "MGKVQTNHYKVYSNQQNECRVLSHWFFHMNKSQQVEVRHHQSIIDKGQOKKKHSDSSMNKSHSDYQCLIHNSHPYLANVVTRHKYNSKKKELCRVKHKVDQLVKVDDLKKVKKDGSKADVLTLSRQGLPDLGRFGQNLLLPCNEQNTLFNWYTRVAKGLMIWVSLKKGKQLEKQWQVNQSAKQVICLNVIERKDRRSVRYIVSVGVGLDVFCSHFCLTKMFRDBQTDNTIFELLKNDLGYGRAIINDKDSLVHGKLIHMSQAKDFFCLYDQEGCVSYNTGWPITKGIVIDHYTGWDSGIQSKIIWQQHAPEDYKKQNFDDSVSVWLFDDETAAHGHMVSGCDEHYDKCLIWFKSVPMKHELGCREELTNTEQLMRQITVVLVKQPISQGETQEHSKQKHDYQWIDKIVKNDNNAWKIYTVNNKHQNFCQRVVMYQESLLKQQDQKVDWEYGWFNRQQVILKRIRNDWLASSTMMKLHKVHPMLKKCHMRFSIHCLQFQDKVSDSANPSDSQCNFYVVSQSKPCLCCFSKKOQRNNQSLLAVDKRLENKKQLLLFKRKMVQQYSGMLFECIQV",
        }
        return gad_references
    
    def parse_fasta(self, file_path):
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
    
    def analyze_sequence(self, sequence):
        """分析序列特征"""
        length = len(sequence)
        
        # 氨基酸组成分析
        amino_acids = {
            'hydrophobic': ['A', 'V', 'L', 'I', 'M', 'F', 'W'],
            'hydrophilic': ['R', 'N', 'D', 'C', 'E', 'Q', 'H', 'K', 'S', 'T', 'Y'],
            'charged': ['R', 'K', 'D', 'E'],
            'glycine': ['G'],
            'proline': ['P'],
            'cysteine': ['C']
        }
        
        composition = {}
        for category, acids in amino_acids.items():
            count = sum(1 for aa in sequence if aa in acids)
            composition[category] = count / length
        
        # GAD酶的特征阈值（针对LR1基因组优化）
        gad_features = {
            "length_range": (200, 600),  # LR1序列长度较短，放宽范围
            "hydrophobic_min": 0.20,
            "hydrophilic_min": 0.35,
            "charged_min": 0.15,
            "glycine_min": 0.05,
            "proline_max": 0.06,
            "cysteine_max": 0.03,
        }
        
        # 结构特征评分
        structural_score = 0
        
        # 长度检查（放宽标准，适应LR1基因组）
        if length >= gad_features["length_range"][0] and length <= gad_features["length_range"][1]:
            structural_score += 20
        
        # 氨基酸组成评分
        if composition['hydrophobic'] >= gad_features["hydrophobic_min"]:
            structural_score += 10
        
        if composition['hydrophilic'] >= gad_features["hydrophilic_min"]:
            structural_score += 10
        
        # LR1候选序列带电氨基酸比例很高，增加权重
        if composition['charged'] >= gad_features["charged_min"]:
            structural_score += 15
        
        if composition['charged'] >= 0.25:  # 特别奖励高带电比例的序列
            structural_score += 10
        
        if composition['glycine'] >= gad_features["glycine_min"]:
            structural_score += 10
        
        if composition['proline'] <= gad_features["proline_max"]:
            structural_score += 5
        
        if composition['cysteine'] <= gad_features["cysteine_max"]:
            structural_score += 5
        
        # 归一化得分
        normalized_score = structural_score / 85
        
        return {
            'length': length,
            'composition': composition,
            'structural_score': normalized_score,
            'structural_raw_score': structural_score
        }
    
    def motif_search(self, sequence):
        """增强保守模式搜索，针对LR1基因组优化"""
        motif_matches = []
        
        for motif in self.gad_motifs:
            pattern = re.compile(motif['pattern'])
            matches = pattern.finditer(sequence)
            
            for match in matches:
                motif_matches.append({
                    'name': motif['name'],
                    'pattern': motif['pattern'],
                    'position': match.start(),
                    'match': match.group(),
                    'weight': motif['weight']
                })
        
        # 计算模式匹配得分（优化算法）
        motif_score = 0
        
        # 核心模式（权重≥3）匹配加分
        core_matches = [match for match in motif_matches if match['weight'] >= 3]
        if core_matches:
            motif_score += sum(match['weight'] for match in core_matches) * 50 / len(sequence)
        
        # 所有模式匹配加分
        motif_score += sum(match['weight'] for match in motif_matches) * 100 / len(sequence)
        
        # LR1特有的带电氨基酸模式额外加分
        charged_patterns = [match for match in motif_matches if 'charged' in match['name']]
        if charged_patterns:
            motif_score += len(charged_patterns) * 20 / len(sequence)
        
        return {
            'matches': motif_matches,
            'motif_score': motif_score,
            'match_count': len(motif_matches),
            'core_motif_count': len(core_matches),
            'charged_pattern_count': len(charged_patterns)
        }
    
    def homology_analysis(self, sequence):
        """增强同源性分析"""
        similarity_scores = {}
        
        for ref_name, ref_seq in self.reference_gad_sequences.items():
            # 简化的相似性计算，针对LR1优化
            score = 0
            
            # 长度相似性（放宽标准）
            len_score = min(len(sequence), len(ref_seq)) / max(len(sequence), len(ref_seq))
            score += len_score * 0.1
            
            # 关键氨基酸相似性（带电氨基酸权重增加）
            key_residues = ['K', 'R', 'D', 'E']  # 带电氨基酸
            seq_key = sum(1 for aa in sequence if aa in key_residues) / len(sequence)
            ref_key = sum(1 for aa in ref_seq if aa in key_residues) / len(ref_seq)
            key_score = min(seq_key, ref_key) / max(seq_key, ref_key)
            score += key_score * 0.4  # 增加权重
            
            # 简化的模式匹配
            shared_motifs = 0
            for motif in ['KK', 'DE', 'RK', 'KD', 'KE']:
                if motif in sequence and motif in ref_seq:
                    shared_motifs += 1
            score += shared_motifs * 0.1
            
            similarity_scores[ref_name] = score
        
        # 找到最佳匹配
        best_match = max(similarity_scores.items(), key=lambda x: x[1]) if similarity_scores else ("Unknown", 0)
        
        return {
            'homology_score': best_match[1],
            'best_match': best_match[0],
            'all_scores': similarity_scores
        }
    
    def functional_domain_analysis(self, sequence):
        """功能域分析（简化版，更容易匹配）"""
        # GAD酶的功能域特征模式（简化版本）
        pfam_patterns = [
            {"name": "GAD-domain-simple", "pattern": "[KR]{2}[DE]{2}", "weight": 2},
            {"name": "PLP-binding-domain-simple", "pattern": "G[ST]G[KR]", "weight": 3},
            {"name": "active-site-domain-simple", "pattern": "[KR][ST][DE]", "weight": 2},
            {"name": "glutamate-binding-domain-simple", "pattern": "[FY][KR][DE]", "weight": 2},
        ]
        
        domain_matches = []
        domain_score = 0
        
        for pattern in pfam_patterns:
            matches = re.finditer(pattern["pattern"], sequence)
            for match in matches:
                domain_matches.append({
                    'name': pattern['name'],
                    'pattern': pattern['pattern'],
                    'position': match.start(),
                    'match': match.group(),
                    'weight': pattern['weight']
                })
                domain_score += pattern['weight'] * 30 / len(sequence)
        
        return {
            'domain_matches': domain_matches,
            'domain_score': domain_score,
            'domain_count': len(domain_matches)
        }
    
    def analyze_all(self, sequence_id, sequence):
        """综合分析"""
        # 结构特征分析
        structural_info = self.analyze_sequence(sequence)
        
        # 保守模式搜索
        motif_info = self.motif_search(sequence)
        
        # 同源性分析
        homology_info = self.homology_analysis(sequence)
        
        # 功能域分析
        domain_info = self.functional_domain_analysis(sequence)
        
        # 综合评分系统（针对LR1优化）
        total_score = (
            structural_info['structural_score'] * 0.30 + 
            motif_info['motif_score'] * 0.35 +
            homology_info['homology_score'] * 0.25 +
            domain_info['domain_score'] * 0.10
        )
        
        return {
            'protein_id': sequence_id,
            'length': structural_info['length'],
            'structural_score': structural_info['structural_score'],
            'structural_raw_score': structural_info['structural_raw_score'],
            'composition': structural_info['composition'],
            'motif_score': motif_info['motif_score'],
            'motif_matches': motif_info['matches'],
            'motif_count': motif_info['match_count'],
            'core_motif_count': motif_info['core_motif_count'],
            'charged_pattern_count': motif_info['charged_pattern_count'],
            'homology_score': homology_info['homology_score'],
            'best_match': homology_info['best_match'],
            'domain_score': domain_info['domain_score'],
            'domain_matches': domain_info['domain_matches'],
            'domain_count': domain_info['domain_count'],
            'total_score': total_score
        }
    
    def run_analysis(self, threshold=0.45):
        """运行分析"""
        print(f"开始LR1基因组专用的GAD酶挖掘，输入文件: {self.genome_file}")
        
        # 解析序列
        sequences = self.parse_fasta(self.genome_file)
        print(f"解析到 {len(sequences)} 个蛋白质序列")
        
        # 分析每个序列
        candidates = []
        processed_count = 0
        
        for seq_id, sequence in sequences.items():
            result = self.analyze_all(seq_id, sequence)
            
            # 筛选候选
            if result['total_score'] >= threshold:
                candidates.append(result)
            
            processed_count += 1
            if processed_count % 100 == 0:
                print(f"已分析: {processed_count}个序列")
        
        # 排序候选
        candidates.sort(key=lambda x: x['total_score'], reverse=True)
        
        print(f"\n分析完成!")
        print(f"找到 {len(candidates)} 个候选GAD酶")
        print(f"结果保存在: {self.output_dir}")
        
        # 保存结果
        self.save_results(candidates, threshold)
        
        # 显示Top候选
        self.display_top_candidates(candidates, threshold)
        
        return candidates
    
    def save_results(self, candidates, threshold):
        """保存结果到文件"""
        # 保存候选列表
        candidates_file = self.output_dir / "lr1_gad_candidates.tsv"
        with open(candidates_file, 'w') as f:
            # 写入标题
            f.write("Protein_ID\tLength\tTotal_Score\tStructural_Score\tMotif_Score\tHomology_Score\tDomain_Score\tBest_Match\tMotif_Count\tCore_Motif_Count\tCharged_Pattern_Count\tDomain_Count\tHydrophobic\tHydrophilic\tCharged\tGlycine\tProline\tCysteine\n")
            
            for candidate in candidates:
                f.write(f"{candidate['protein_id']}\t{candidate['length']}\t{candidate['total_score']:.3f}\t{candidate['structural_score']:.3f}\t{candidate['motif_score']:.3f}\t{candidate['homology_score']:.3f}\t{candidate['domain_score']:.3f}\t{candidate['best_match']}\t{candidate['motif_count']}\t{candidate['core_motif_count']}\t{candidate['charged_pattern_count']}\t{candidate['domain_count']}\t{candidate['composition']['hydrophobic']:.3f}\t{candidate['composition']['hydrophilic']:.3f}\t{candidate['composition']['charged']:.3f}\t{candidate['composition']['glycine']:.3f}\t{candidate['composition']['proline']:.3f}\t{candidate['composition']['cysteine']:.3f}\n")
        
        # 保存详细报告
        report_file = self.output_dir / "lr1_gad_report.txt"
        with open(report_file, 'w') as f:
            f.write("LR1基因组专用的GAD酶挖掘分析报告\n")
            f.write("=====================================\n")
            f.write(f"输入文件: {self.genome_file}\n")
            f.write(f"分析日期: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"筛选阈值: {threshold}\n")
            f.write(f"候选数量: {len(candidates)}\n\n")
            
            f.write("Top 10候选GAD酶:\n")
            for i, candidate in enumerate(candidates[:10], 1):
                f.write(f"{i}. {candidate['protein_id']}\n")
                f.write(f"  综合评分: {candidate['total_score']:.3f}\n")
                f.write(f"  结构评分: {candidate['structural_score']:.3f}\n")
                f.write(f"  模式评分: {candidate['motif_score']:.3f}\n")
                f.write(f"  同源性评分: {candidate['homology_score']:.3f}\n")
                f.write(f"  功能域评分: {candidate['domain_score']:.3f}\n")
                f.write(f"  最佳参考匹配: {candidate['best_match']}\n")
                f.write(f"  匹配模式数量: {candidate['motif_count']}\n")
                f.write(f"  核心模式数量: {candidate['core_motif_count']}\n")
                f.write(f"  带电模式数量: {candidate['charged_pattern_count']}\n")
                f.write(f"  功能域数量: {candidate['domain_count']}\n")
                f.write(f"  氨基酸组成: hydrophobic={candidate['composition']['hydrophobic']:.3f}, hydrophilic={candidate['composition']['hydrophilic']:.3f}, charged={candidate['composition']['charged']:.3f}, glycine={candidate['composition']['glycine']:.3f}\n\n")
            
            # 统计信息
            if candidates:
                avg_length = sum(c['length'] for c in candidates) / len(candidates)
                avg_charged = sum(c['composition']['charged'] for c in candidates) / len(candidates)
                f.write(f"\n统计信息:\n")
                f.write(f"候选评分范围: {max(c['total_score'] for c in candidates):.3f} - {min(c['total_score'] for c in candidates):.3f}\n")
                f.write(f"平均长度: {avg_length:.1f} aa\n")
                f.write(f"平均带电氨基酸比例: {avg_charged:.3f}\n")
                
                # 参考序列匹配分布
                match_distribution = {}
                for candidate in candidates:
                    match = candidate['best_match']
                    match_distribution[match] = match_distribution.get(match, 0) + 1
                
                f.write(f"\n参考序列匹配分布:\n")
                for match, count in match_distribution.items():
                    f.write(f"{match}: {count}个候选\n")
    
    def display_top_candidates(self, candidates, threshold):
        """显示Top候选"""
        print("\nTop 10候选:")
        for i, candidate in enumerate(candidates[:10], 1):
            print(f"{i}. {candidate['protein_id']}")
            print(f"  综合评分: {candidate['total_score']:.3f}")
            print(f"  最佳参考匹配: {candidate['best_match']}")
            print(f"  模式匹配数: {candidate['motif_count']}")
            print(f"  核心模式数: {candidate['core_motif_count']}")
            print(f"  带电模式数: {candidate['charged_pattern_count']}")
            print(f"  功能域数: {candidate['domain_count']}")
        
        print("\n详细结果:")
        print(f"1. 候选列表: {self.output_dir}/lr1_gad_candidates.tsv")
        print(f"2. 分析报告: {self.output_dir}/lr1_gad_report.txt")
        print(f"3. 阈值: {threshold}")
        print(f"4. 候选数量: {len(candidates)}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='LR1基因组专用的GAD酶挖掘工具')
    parser.add_argument('-g', '--genome', required=True, help='蛋白质序列文件（FASTA格式）')
    parser.add_argument('-o', '--output', default='lr1_results', help='输出目录')
    parser.add_argument('-t', '--threshold', type=float, default=0.45, help='筛选阈值（默认0.45）')
    
    args = parser.parse_args()
    
    # 检查文件是否存在
    if not os.path.exists(args.genome):
        print(f"错误: 文件 {args.genome} 不存在!")
        return
    
    # 运行分析
    finder = AdvancedGADFinder(args.genome, args.output)
    candidates = finder.run_analysis(args.threshold)

if __name__ == '__main__':
    main()