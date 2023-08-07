import pandas


class_version = "V2.1(Editor) 2023-08-07"


# Class
# 声明一个类，存储所有gene的信息
class Total(object):
    def __init__(self, sample_list, gene_dict={}):
        self.gene_dict = gene_dict  # 存储原始的Gene对象
        self.sample_list = sample_list
        self.start_dict = {}  # 为了加快check_gene_id()的检索速度测试使用 {<start>: [<gene_id>, <gene_id>],
                              #                                          <start>: [<gene_id>, <gene_id>, <gene_id>, ...]}

    def __check_gene_id(self, Gene):
        # change:
        #   检验gene是否已被记录
        #       1.gene_id是否已被记录
        #       2.在相应染色体及正负链上, 该gene的范围是否已被记录
        # output:
        #   str, 若该gene_id已被记录则返回'0', 若该gene_id未被记录(start未被记录)则返回'1',
        #        若该gene_id未被记录(start已被记录但end未被记录)则返回'2'
        #        若该gene_id未被记录(start与end均已被记录)则返回被记录的gene_id

        if Gene.gene_id in self.gene_dict.keys():
            return '0'
        else:
            if Gene.start not in self.start_dict.keys():
                return '1'
            else:
                gene_id_list = self.start_dict.get(Gene.start)
                for existed_gene in gene_id_list:
                    existed_gene = self.gene_dict.get(existed_gene)
                    if existed_gene.strand != Gene.strand:
                        continue
                    elif existed_gene.chr != Gene.chr:
                        continue
                    else:
                        if existed_gene.end == Gene.end:
                            return existed_gene.gene_id
                return '2'

    def add_gene(self, Gene):
        # change:
        #   若gene信息未被记录, 则向gene_dict添加Gene对象，并向start_dict添加相关信息
        # output:
        #   str, 若成功添加该gene信息则返回"True"，若该gene_id已被记录则返回"False",
        #        若该gene_id未被记录但相关信息已被记录则返回相应的gene_id

        check_exist = self.__check_gene_id(Gene)
        if check_exist == '0':
            return "False"
        elif check_exist == '1':
            self.gene_dict[Gene.gene_id] = Gene
            self.start_dict[Gene.start] = [Gene.gene_id]
            return "True"
        elif check_exist == '2':
            self.gene_dict[Gene.gene_id] = Gene
            temp = self.start_dict.get(Gene.start)
            temp.append(Gene.gene_id)
            self.start_dict[Gene.start] = temp
            return "True"
        else:
            return check_exist

    def get_df(self, sample_name_list=None):
        '''
        input:
            sample_name_list,\t list,\t 输出的df中所要包含的统计样本\n
        change:
            整理该对象下所有gene的信息到一个df中
        output:
            df,\t pandas.DataFrame,\t 含有列chr, strand, source, gene_id, gene_name, gene_biotype,
                                        transcript_id, transcript_name, transcript_biotype,
                                        transcript_start, transcript_end, ...\n
        '''
        sampleNameList = (sample_name_list, [])[sample_name_list is None]

        df = {}
        columns = ["chr", "strand", "source", "gene_id", "gene_name", "gene_biotype",
                "transcript_id", "transcript_name", "transcript_biotype", "transcript_start", "transcript_end"]
        columns = columns + ["{}_counts".format(sample) for sample in sampleNameList]
        columns = columns + ["{}_relative".format(sample) for sample in sampleNameList]

        # 获取所有gene中所有transcript的信息
        for geneId, geneObject in self.gene_dict.items():
            df.update(geneObject.getTranscriptInfo(sampleNameList=sample_name_list))
        # 将dict转为df
        df = pandas.DataFrame.from_dict(df, orient="index", columns=columns)
        # 美化格式
        df = df.sort_values(by=["strand", "chr"])

        return df

    def get_exon_combination(self):
        """
        input:
        change:
            统计Total对象中每一个gene的exon
        output:
            exon_combination, dict, {"<gene_id1>": ["<exon1_start>-<exon1_end>",
                                                    "<exon2_start>-<exon2_end>", ...
                                                    ],
                                     "<gene_id2>": ["<exon1_start>-<exon1_end>",
                                                    "<exon2_start>-<exon2_end>", ...
                                                    ], ...
                                     } 
        """
        exon_combination = {}
        for gene_id in self.gene_dict.keys():
            gene_object = self.gene_dict.get(gene_id)
            temp_exon_combination = gene_object.get_exon_combination()
            exon_combination.update(temp_exon_combination)

        return exon_combination

    def add_gene_classified(self):
        """
        change:
            遍历每一个gene对象，使用对象中的add_gene_classified方法
                向对象添加gene_classified属性
        output:
            dict,\t {<gene_id>: {"range_start": set(range_start1, range_start2, ...),
                                 "range_end": set(range_end1, range_end2, ...)}
                     <gene_id>: ...}
        """
        TSS_TES_dict = {}
        for gene_id in self.gene_dict.keys():
            TSS_TES_dict[gene_id] = self.gene_dict[gene_id].add_gene_classified()

        return TSS_TES_dict

    def computeRelativeExpression(self):
        sampleList = self.sample_list  # list, 存储所有样本的样本名

        # 统计每个样本中所有transcript的rep值的和, 以供后续计算样本中的相对表达值
        countsSample = {sample: 0 for sample in sampleList}
        for geneId, geneObject in self.gene_dict.items():
            # 遍历每一个gene
            for transcriptId, transcriptObject in geneObject.transcript_dict.items():
                # 遍历每一个transcript
                for sample, rep in transcriptObject["countsExpression"].items():
                    # 遍历transcript在每个样本中的计数
                    countsSample[sample] = countsSample.get(sample, 0) + rep

        # 对每个transcript计算相对表达值
        for geneId, geneObject in self.gene_dict.items():
            # 遍历每一个gene
            for transcriptId, transcriptObject in geneObject.transcript_dict.items():
                # 遍历每一个transcript
                for sample, rep in transcriptObject["countsExpression"].items():
                    # 遍历transcript在每个样本中的计数
                    # 计算每个transcript在每个样本中的相对表达值
                    # 单个transcript在单个样本中的相对表达值 = 10^6 * 单个transcript在单个样本中的计数/单个样本中所有transcript的计数的和
                    self.gene_dict[geneId].transcript_dict[transcriptId]["relativeExpression"][sample] = rep / countsSample[sample] * 10**6

        return None


# 定义一个类，以存储基因的相关信息
#   具有属性：染色体号，source，正负链，start/end，包含的转录本id
#   存储：转录本名称及范围、外显子范围
#   具有方法：
#       1.添加外显子
#       2.添加转录本
#       3.统计该基因的转录本数量

# 需要建立一个df, 可以实现gene_id与transcript_id的互查

# 类
class Gene(object):
    def __init__(self, chr, gene_id, gene_name, gene_biotype,
                 source, strand, start, end,
                 transcript_dict, exon_dict):
        self.chr = chr  # str
        self.gene_id = gene_id  # str
        self.gene_name = gene_name  # str
        self.gene_biotype = gene_biotype  # str  ["protein", "non_protein", "un_classified"]
        self.source = source  # str
        self.strand = strand  # str
        self.start = int(start)  # int
        self.end = int(end)  # int
        self.transcript_dict = transcript_dict  # {<transcript_id>: {"transcript_name": <transcript_name>
                                                #                    "transcript_biotype": <transcript_biotype="un_classified">,
                                                #                    "range": [<start>, <end>],
                                                #                    "exon_range": {<start>: [<end>], <start>: [<end>, <end>], ...},
                                                #                    "countsExpression" : {<sample_name1>: <expression>,
                                                #                                          <sample_name2>: <expression>}
                                                #                    ！"relativeExpression": {<sample_name1>: <expression>,
                                                #                                           <sample_name2>: <expression>, ...}
                                                #  <transcript_id>: ...,
                                                #  ...
                                                # }
        self.exon_dict = exon_dict  # {start: [end], start: [end, end, ...], ...}  int
        # 初始化属性
        self.gene_classified = None

    def __check_exon_exist(self, exon_start, exon_end):
        # change:
        #   检查指定的exon是否已被记录
        # output:
        #   int, 若start未被记录则返回0, 若start被记录end未被记录则返回1, 若该exon的start与end均已被记录则返回2

        if exon_start not in self.exon_dict.keys():
            return 0
        else:
            if exon_end not in self.exon_dict.get(exon_start):
                return 1
            else:
                return 2

    def add_exon(self, exon_start, exon_end):
        # change:
        #   对基因增加一个exon
        # output:
        #   str, 若成功添加exon则返回"success_add", 若该exon已被记录则返回"existed_exon"

        exon_exist_mark = self.__check_exon_exist(exon_start, exon_end)
        if exon_exist_mark == 0:
            # 存储start与end未记录的exon的信息
            self.exon_dict[exon_start] = [exon_end]
            return "success add"
        elif exon_exist_mark == 1:
            # 存储start已记录而end未记录的exon的信息
            self.exon_dict[exon_start].append(exon_end)
            return "success add"
        else:
            return "existed_exon"

    def __check_transcript_exist(self, transcript_id, start, end, exon_list):
        # change:
        #   检验指定的transcript是否已被记录
        # output:
        #   str, 若trnascript_id已被记录则返回transcript_id, 表示该transcript已被记录, 若transcript_id未被记录则检验该transcript的range是否已被记录
        #           若该transcript的range未被记录，则返回1，表示该transcript未被记录
        #           若该transcript的range已被记录, 则进一步检验该transcript的exon组成是否已被记录
        #               若该transcript的exon组成已被记录, 则返回exist_transcript_id, 表示该transcript已被记录
        #               若该transcript的exon组成未被记录，则返回3, 表示该transcript未被记录

        if transcript_id in self.transcript_dict.keys():
            return transcript_id  # transcript_id已被记录
        else:
            for exist_transcript_id in self.transcript_dict.keys():
                if start != self.transcript_dict[exist_transcript_id]["range"][0]:
                    continue
                else:
                    if end != self.transcript_dict[exist_transcript_id]["range"][1]:
                        continue
                    else:
                        # transcript的id未被记录, 但start与end均已被记录且与exist_transcript_id的start与end一致
                        # 接下来比较transcript与exist_transcript的exon组成是否一致

                        # 将exist_transcript的exon_range由dict转为list
                        exist_transcript_exon = []
                        for exist_start,exist_end_list in self.transcript_dict[exist_transcript_id]["exon_range"].items():
                            for exist_end in exist_end_list:
                                exist_transcript_exon.append([exist_start, exist_end])

                        # 判断exon的数量是否一致,
                        length = len(exon_list)
                        if length != len(exist_transcript_exon):
                            # 若不一致则表示该transcript未被记录
                            return '3'  # range一致但exon组成不一致, 新transcript
                        else:
                            # 若一致，判断exon组成是否一致
                            # 将两个exon_list按start的位置升序排序，在两两比较list是否一致
                            exon_list = sorted(exon_list, key=lambda x: x[0])
                            exist_transcript_exon = sorted(exist_transcript_exon, key=lambda x: x[0])
                            for i in range(0,length):
                                if exon_list[i] != exist_transcript_exon[i]:
                                    # 若存在不一致，则表示transcript为新transcript
                                    return '3'  # range一致但exon组成不一致, 新transcript
                            return exist_transcript_id  # range一致且exon一致, 旧transcript
            return '1'  # range未被记录, 新transcript

    def add_transcript(self, transcript_id, transcript_name, transcript_biotype,
                       start, end, exon_list, sample_name, sample_counts):
        # change:
        #   对基因增加一个新的transcript, 并记录该transcript在相应样本中的counts数
        # output:
        #   bool, 若成功添加该transcript并添加了相应counts数则返回True，若未添加transcript只添加了相应counts数则返回False

        exist_mark = self.__check_transcript_exist(transcript_id, start, end, exon_list)
        if exist_mark == '1' or exist_mark == '3':
            # 将exon_list转换为transcript_dict中exon_range的格式  {<start>: [<end>], <start>: [<end>, <end>], ...}
            exon_range = {}
            for [exon_start, exon_end] in exon_list:
                if exon_start not in exon_range.keys():
                    exon_range[exon_start] = [exon_end]  # list type
                else:
                    temp = exon_range[exon_start]
                    temp.append(exon_end)
                    exon_range[exon_start] = temp
            # gene中的transcript_dict新增一个transcript_id，并添加相应的键值对
            self.transcript_dict[transcript_id] = {"transcript_name": transcript_name,
                                                   "transcript_biotype": transcript_biotype,
                                                   "range": [start, end],
                                                   "exon_range": exon_range,
                                                   "countsExpression": {sample_name: sample_counts},
                                                   "relativeExpression": {},
                                                   }
            return True
        else:
            # 仅更新exist_transcript_id在相应样本中的counts数
            self.transcript_dict[exist_mark]["countsExpression"][sample_name] = sample_counts
            return False

    def getTranscriptInfo(self, geneInfo=True, sampleNameList=None):
        '''
        input:
            geneInfo,\t bool,\t, 是否要在dict中保留gene的相关信息(strand, gene_id等)
            sampleNameList,\t list,\t 需要进行展示的统计样本的list
        change:
            整理该gene中的所有信息，返回一个dict
        output:
            info, dict, {<transcriptId>: {"chr": <chr>, "strand": <strand>, "source": <source>,
                                        "gene_id": <gene_id>, "gene_name": <gene_name>, "gene_biotype": <gene_biotype>,
                                        "transcript_id": <transcript_id>, "transcript_name": <transcript_name>, "transcript_biotype", <transcript_biotype>,
                                        ...},
                        <transcriptId>: ...}
        '''
        geneInfo = geneInfo
        sampleNameList = (sampleNameList, [])[sampleNameList is None]

        info = {}
        for transcriptId, transcriptObject in self.transcript_dict.items():
            info[transcriptId] = {"transcript_id": transcriptId,
                                "transcript_name": transcriptObject["transcript_name"],
                                "transcript_biotype": transcriptObject["transcript_biotype"],
                                "transcript_start": transcriptObject["range"][0],
                                "transcript_end": transcriptObject["range"][1],
                                }
            if geneInfo is True:
                info[transcriptId].update({"chr": self.chr,
                                        "strand": self.strand,
                                        "source": self.source,
                                        "gene_id": self.gene_id,
                                        "gene_name": self.gene_name,
                                        "gene_biotype": self.gene_biotype,
                                        })
            for sample_name in sampleNameList:
                info[transcriptId]["{}_counts".format(sample_name)] = transcriptObject["countsExpression"].get(sample_name, 0)
                info[transcriptId]["{}_relative".format(sample_name)] = transcriptObject["relativeExpression"].get(sample_name, 0)

        return info

    def get_exon_combination(self):
        """
        input:
        change:
            统计该gene中的全部exon。遍历self.exon_dict, 得到exon, 以"<start>-<end>"形式返回每一个exon的位置.
        output:
            exon_combination, dict, {"<gene_id>": ["<exon1_start>-<exon1_end>",
                                                   "<exon2_start>-<exon2_end>", ...
                                                   ]
                                     }
        """
        exon_combination = {}
        for exon_start, exon_end_list in self.exon_dict.items():
            for exon_end in exon_end_list:
                temp = exon_combination.get(self.gene_id, [])
                temp.append("{}-{}".format(exon_start, exon_end))
                exon_combination[self.gene_id] = temp

        return exon_combination

    def get_transcript_according_exon(self, exon_checked):
        """
        input:
            exon_checked, str, "<exon_start>-<exon_end>"
        change:
            查询gene中具有指定exon的transcript
        output:
            transcript_checked, list, ["transcript_id", ...]
        """
        exon_checked = exon_checked

        exon_checked = [int(i) for i in exon_checked.split('-')]
        exon_check_start = exon_checked[0]
        exon_check_end = exon_checked[1]

        transcript_checked = []
        for transcript_id, transcript_info_dict in self.transcript_dict.items():
            if exon_check_start in transcript_info_dict["exon_range"].keys():
                end_list = transcript_info_dict["exon_range"][exon_check_start]
                if exon_check_end in end_list:
                    transcript_checked.append(transcript_id)

        return transcript_checked

    def add_gene_classified(self):
        """
        change:
            根据gene中所具有的转录起始位点及转录终止位点的数量, 确定gene的类型(TSS-PAS, TSS-APA, ATSS-PAS, ATSS-APA)
                将其类型作为gene的属性gene_classified添加到对象中
            注意：
                1.已考虑正负链
        output:
            dict,\t {"range_start": range_start, "range_end": range_end}
        """

        range_start = set()
        range_end = set()
        for vdict in self.transcript_dict.values():
            transcript_range = vdict.get("range")
            # 根据gene的strand, 确定gene中存在的转录起始位点以及转录终止位点的位置
            if self.strand == '-':
                range_start.add(transcript_range[1])
                range_end.add(transcript_range[0])
            else:
                range_start.add(transcript_range[0])
                range_end.add(transcript_range[1])

        # 根据gene中所具有的转录起始位点及转录终止位点的数量，确定gene的类型(TSS-PAS, TSS-APA, ATSS-PAS, ATSS-APA)
        self.gene_classified = "{}-{}".format(("ATSS", "TSS")[len(range_start)==1],
                                              ("APA", "PAS")[len(range_end)==1]
                                              )
        '''
        if len(range_start) == 1:
            if len(range_end) == 1:
                self.gene_classified = "TSS-PAS"
            else:
                self.gene_classified = "TSS-APA"
        else:
            if len(range_end) == 1:
                self.gene_classified = "ATSS-PAS"
            else:
                self.gene_classified = "ATSS-APA"
        '''

        return {"range_start": range_start, "range_end": range_end}
