import pandas


class_version = "V1.0(Editor) 2023-07-13"


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


    def get_df(self, sample_name_list=[]):
        # change:
        #   整理该对象下所有gene的信息到一个df中
        # output:
        #   df, pandas.DataFrame, 含有列chr, strand, source, gene_id, gene_name, gene_biotype,
        #                               transcript_id, transcript_name, transcript_biotype, transcript_start, transcript_end, sample_name...
        sample_name_list = sample_name_list

        gene_id_list = list(self.gene_dict.keys())
        columns = ["chr", "strand", "source",
                   "gene_id", "gene_name", "gene_biotype",
                   "transcript_id", "transcript_name", "transcript_biotype",
                   "transcript_start", "transcript_end"]
        columns = columns + sample_name_list
        df = pandas.DataFrame(columns=columns)

        # 为了加快concat速度
        remain_num = len(gene_id_list)  # 进度条
        while remain_num != 0:
            print("{}".format(remain_num), end='\r')  # 打印进度条

            if remain_num >= 30:
                temp1 = gene_id_list.pop()
                temp2 = gene_id_list.pop()
                temp3 = gene_id_list.pop()
                temp4 = gene_id_list.pop()
                temp5 = gene_id_list.pop()
                temp6 = gene_id_list.pop()
                temp7 = gene_id_list.pop()
                temp8 = gene_id_list.pop()
                temp9 = gene_id_list.pop()
                temp10 = gene_id_list.pop()
                temp11 = gene_id_list.pop()
                temp12 = gene_id_list.pop()
                temp13 = gene_id_list.pop()
                temp14 = gene_id_list.pop()
                temp15 = gene_id_list.pop()
                temp16 = gene_id_list.pop()
                temp17 = gene_id_list.pop()
                temp18 = gene_id_list.pop()
                temp19 = gene_id_list.pop()
                temp20 = gene_id_list.pop()
                temp21 = gene_id_list.pop()
                temp22 = gene_id_list.pop()
                temp23 = gene_id_list.pop()
                temp24 = gene_id_list.pop()
                temp25 = gene_id_list.pop()
                temp26 = gene_id_list.pop()
                temp27 = gene_id_list.pop()
                temp28 = gene_id_list.pop()
                temp29 = gene_id_list.pop()
                temp30 = gene_id_list.pop()

                temp1 = self.gene_dict[temp1].get_df(sample_name_list)
                temp2 = self.gene_dict[temp2].get_df(sample_name_list)
                temp3 = self.gene_dict[temp3].get_df(sample_name_list)
                temp4 = self.gene_dict[temp4].get_df(sample_name_list)
                temp5 = self.gene_dict[temp5].get_df(sample_name_list)
                temp6 = self.gene_dict[temp6].get_df(sample_name_list)
                temp7 = self.gene_dict[temp7].get_df(sample_name_list)
                temp8 = self.gene_dict[temp8].get_df(sample_name_list)
                temp9 = self.gene_dict[temp9].get_df(sample_name_list)
                temp10 = self.gene_dict[temp10].get_df(sample_name_list)
                temp11 = self.gene_dict[temp11].get_df(sample_name_list)
                temp12 = self.gene_dict[temp12].get_df(sample_name_list)
                temp13 = self.gene_dict[temp13].get_df(sample_name_list)
                temp14 = self.gene_dict[temp14].get_df(sample_name_list)
                temp15 = self.gene_dict[temp15].get_df(sample_name_list)
                temp16 = self.gene_dict[temp16].get_df(sample_name_list)
                temp17 = self.gene_dict[temp17].get_df(sample_name_list)
                temp18 = self.gene_dict[temp18].get_df(sample_name_list)
                temp19 = self.gene_dict[temp19].get_df(sample_name_list)
                temp20 = self.gene_dict[temp20].get_df(sample_name_list)
                temp21 = self.gene_dict[temp21].get_df(sample_name_list)
                temp22 = self.gene_dict[temp22].get_df(sample_name_list)
                temp23 = self.gene_dict[temp23].get_df(sample_name_list)
                temp24 = self.gene_dict[temp24].get_df(sample_name_list)
                temp25 = self.gene_dict[temp25].get_df(sample_name_list)
                temp26 = self.gene_dict[temp26].get_df(sample_name_list)
                temp27 = self.gene_dict[temp27].get_df(sample_name_list)
                temp28 = self.gene_dict[temp28].get_df(sample_name_list)
                temp29 = self.gene_dict[temp29].get_df(sample_name_list)
                temp30 = self.gene_dict[temp30].get_df(sample_name_list)
                df = pandas.concat([df, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10,
                                    temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20,
                                    temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28, temp29, temp30],
                                    axis=0)
            else:
                temp1 = gene_id_list.pop()
                temp1 = self.gene_dict[temp1].get_df(sample_name_list)
                df = pandas.concat([df, temp1], axis=0)

            remain_num = len(gene_id_list)

        return df


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
        # 可删除 self.counts_dict = counts_dict  # {sample_name1: counts, sample_name2: counts, ...}
        self.transcript_dict = transcript_dict  # {<transcript_id>: {"transcript_name": <transcript_name>
                                                #                    "transcript_biotype": <transcript_biotype="un_classified">,
                                                #                    "range": [<start>, <end>],
                                                #                    "exon_range": {<start>: [<end>], <start>: [<end>, <end>], ...},
                                                #                    <sample_name1>: <counts>,
                                                #                    <sample_name2>: <counts>, ...},
                                                #  <transcript_id>: ...,
                                                #  ...
                                                # }
        self.exon_dict = exon_dict  # {start: [end], start: [end, end, ...], ...}

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
                                                   sample_name: sample_counts
                                                   }
            return True
        else:
            # 仅更新exist_transcript_id在相应样本中的counts数
            self.transcript_dict[exist_mark][sample_name] = sample_counts
            return False


    def get_df(self, sample_name_list=[]):
        # change:
        #   整理该gene中的所有信息，返回一个pandas.DataFrame对象
        # output:
        #   df, pandas.DataFrame, 含有列chr, strand, source, gene_id, transcript_id, transcript_start, transcript_end, sample_name...
        sample_name_list = sample_name_list

        chr = self.chr
        strand = self.strand
        source = self.source
        gene_id = self.gene_id
        gene_name = self.gene_name
        gene_biotype = self.gene_biotype

        columns = ["chr", "strand", "source",
                   "gene_id", "gene_name", "gene_biotype",
                   "transcript_id", "transcript_name", "transcript_biotype",
                   "transcript_start", "transcript_end"]
        columns = columns + sample_name_list
        df = pandas.DataFrame(index=list(self.transcript_dict.keys()),
                              columns=columns)

        for transcript_id in self.transcript_dict.keys():
            transcript_name = self.transcript_dict[transcript_id]["transcript_name"]
            transcript_biotype = self.transcript_dict[transcript_id]["transcript_biotype"]
            transcript_start = self.transcript_dict[transcript_id]["range"][0]
            transcript_end = self.transcript_dict[transcript_id]["range"][1]

            df.at[transcript_id, "chr"] = chr
            df.at[transcript_id, "strand"] = strand
            df.at[transcript_id, "source"] = source
            df.at[transcript_id, "gene_id"] = gene_id
            df.at[transcript_id, "gene_name"] = gene_name
            df.at[transcript_id, "gene_biotype"] = gene_biotype
            df.at[transcript_id, "transcript_id"] = transcript_id
            df.at[transcript_id, "transcript_name"] = transcript_name
            df.at[transcript_id, "transcript_biotype"] = transcript_biotype
            df.at[transcript_id, "transcript_start"] = transcript_start
            df.at[transcript_id, "transcript_end"] = transcript_end

            for sample_name in sample_name_list:
                df.at[transcript_id, sample_name] = self.transcript_dict[transcript_id].get(sample_name, 0)

        return df
