import pandas
import time
import re
import os
import gzip
from def_class import *


__version__ = "V2.5(Editor) 2023-08-22"


# Function
def log(function):
    name = function.__name__
    def wrapper(*args, **kwargs):
        tab = kwargs.get("tab_level", 0)

        [print('\n') if tab == 0 else None]

        print("{}[Function]{} start.".format('\t'*tab, name))
        print("{}[Time]{}".format('\t'*(tab+1),
                                  time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))))
        for key, value in kwargs.items():
            if type(value) in [int, float, str, bool, list]:
                print("{}[Paraments]{}: {}".format('\t'*(tab+1),
                                                   key, value))
            else:
                print("{}[Paraments]{}: <...>".format('\t'*(tab+1),
                                                      key))
        result = function(*args, **kwargs)
        print("{}[{}]{} finished.".format('\t'*tab,
                                          time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                          name))
        return result
    return wrapper


def get_file_path(file_path, sample_name_list, tab_level=0):
    """
    input:
        file_path, str, data文件夹路径
        sample_name_list, list, 含有样本名的列表
        tab_level, int, the parament for log format, there is no need to change
    change:
        搜索file_path下的每个样本文件夹中的文件名,
            记录annotation及quantification文件的绝对路径
    output:
        sample_file_location, dict, 该字典存储每个样本的annotation及quantification文件的绝对路径
                                    {<sample1>:{"annotation": <dir1>, "quantification": <dir2>},
                                     <sample2>:{"annotation": <dir1>, "quantification": <dir2>}, ...}
    """

    file_path = file_path
    sample_name_list = sample_name_list

    sample_file_location = {sample: {} for sample in sample_name_list}
    for sample in sample_name_list:
        filename_list = os.listdir("{}raw_data/{}".format(file_path, sample))
        for filename in filename_list:
            if re.search(pattern="annotation", string=filename) is not None:
                sample_file_location[sample]["annotation"] = "{}raw_data/{}/{}".format(file_path, sample, filename)
            elif re.search(pattern="quantification", string=filename) is not None:
                sample_file_location[sample]["quantification"] = "{}raw_data/{}/{}".format(file_path, sample, filename)
            else:
                continue

    return sample_file_location


def load_annotation(input_filename, additional_info=["gene_id", "transcript_id", "gene_name"],
                    only_type="transcript", tab_level=0):
    input_filename = input_filename
    additional_info_list = additional_info
    only_type = only_type

    columns_name = ["chr", "source", "type", "start", "end", "strand"]
    columns_name = columns_name + additional_info_list

    df = []
    with open(input_filename, 'r') as file:
        for line in file:
            line = line.strip('\n')
            line = line.split('\t')

            # 只保留transcript
            if only_type == "all":
                pass
            elif line[2] != only_type:
                continue

            # 过滤掉 chr为ERCC 或 未能准确匹配到染色体上 的注释
            if line[0][0:4] == "ERCC":
                continue
            elif line[0][0:5] == "chrUn":
                continue
            elif line[0][-6:] == "random":
                continue

            # 过滤掉start==end的行(gene, transcript, exon)
            if line[3] == line[4]:
                continue

            # 去除不想要的列
            line.pop(7)
            line.pop(5)

            # 统计最后一列的信息
            additional_info = {}
            temp_line = line[-1].split(';')
            temp_line.remove('')
            temp_line = [i.strip(' ') for i in temp_line]
            for i in temp_line:
                key,value = i.split(' ')
                value = value.strip('"')
                additional_info[key] = value

            # 保存前8列以及指定信息
            line = line[0:6]
            for i in additional_info_list:
                line.append(additional_info.get(i,""))

            df.append(line)

    df = pandas.DataFrame(df, columns=columns_name)

    return df


def load_quantification(input_filename,
                        sample_name,
                        cutoff_value=1,
                        col=["annot_gene_id", "annot_transcript_id", "gene_novelty", "transcript_novelty"],
                        tab_level=0):
    """
    input:
        input_filename, str, the absolute path of quantification file
        sample_name, str, the name of sample
        col, list, the columns need to be retained of the df
    change:
        load input file and retain the specified columns
        filter rows according to cutoff_value
    output:
    """
    input_filename = input_filename
    sample_name = sample_name
    cutoff_value = cutoff_value
    col = col

    df = pandas.read_csv(input_filename, sep='\t')
    col_list = df.columns.to_list()
    df.rename(columns={col_list[-1]: sample_name}, inplace=True)
    col = [col_list.index(x) for x in col]
    col.append(len(col_list)-1)  # the last column in df
    df = df.iloc[:, col]

    sig_row = df[sample_name].map(lambda x: x >= cutoff_value)
    df = df.loc[sig_row, :]

    return df


def merge_quantification_annotation(quantification, annotation,
                                    tab_level=0):
    df_quan = quantification
    df_annotation = annotation

    # 注意：这一步删除了annotation中不存在于quantification的transcript
    df = pandas.merge(df_quan, df_annotation,
                      left_on="annot_transcript_id", right_on="transcript_id",
                      how="left")
    # 删除quantification中不存在于annotation的transcript
    df = df.dropna(axis=0, how="any", subset="transcript_id")
    # 删除df中的gene_id列与transcript_id列，因为与quantification重复
    df = df.drop(labels=["gene_id", "transcript_id"], axis=1)

    return df


def get_gene_info(df, tab_level=0):
    # 统计每个基因的TSS及TES的数量
    df = df

    df_group = df.groupby("annot_gene_id")
    gene_info = pandas.DataFrame(index=list(set(df["annot_gene_id"].to_list())),
                                columns=["start","end","strand","transcripts"])

    process_n = 0
    process_end = len(df_group)
    for gene_id in df_group.groups.keys():
        # 进度条
        process_n += 1
        print("\t[Process]{}/{}".format(process_n, process_end), end='\r')
        # 获得单个基因的所有转录本信息
        temp_df = df_group.get_group(gene_id)

        num_transcript = temp_df.shape[0]
        strand = list(set(temp_df["strand"].to_list()))
        # 根据strad对起始点及终止点的数量进行统计
        if len(strand) == 1 and strand[0] == '-':
            # 负链
            num_end = len(set(temp_df["start"].to_list()))
            num_start = len(set(temp_df["end"].to_list()))
            strand = strand[0]
        elif len(strand) == 1 and strand[0] == '+':
            # 正链
            num_start = len(set(temp_df["start"].to_list()))
            num_end = len(set(temp_df["end"].to_list()))
            strand = strand[0]
        else:
            raise KeyError("[错误]strand: {}".format(strand))

        gene_info.at[gene_id, "start"] = num_start
        gene_info.at[gene_id, "end"] = num_end
        gene_info.at[gene_id, "strand"] = strand
        gene_info.at[gene_id, "transcripts"] = num_transcript

    return gene_info


@log
def load_sample(total, annotation_location, quantification_location, sample_name, 
                counts_cutoff, annotation_col, quantification_col,
                tab_level=0):
    """
    input:
        total, Total, 含有整理的所有信息的对象
        annotation_location, str, 单个样本的annotation文件
        quantification_location, str, 单个样本的quantification文件
        sample_name, str, 单个样本的样本名称
        counts_cutoff, int, quantification文件中counts的阈值(>=value)
        annotation_col, list, 在annotation文件中要保留的额外的列名
        quantification_col, list, 在quantification文件中要保留的额外的列名
    change:
        整合单个样本的annotation和quantification数据到total对象中
    output:
        total, Total, 存储整合annotation和quantification文件得到的所有数据
    """

    total = total
    annotation_location = annotation_location
    quantification_location = quantification_location
    sample_name = sample_name
    counts_cutoff = counts_cutoff
    annotation_col = annotation_col
    quantification_col = quantification_col

    # 获取单个样本的transcript信息
    df_annotation = load_annotation(input_filename=annotation_location,
                                    additional_info=annotation_col,
                                    only_type="transcript",
                                    tab_level=tab_level+1)
    df_quantification = load_quantification(input_filename=quantification_location,
                                            sample_name=sample_name,
                                            cutoff_value=counts_cutoff,
                                            col=quantification_col,
                                            tab_level=tab_level+1)
    df_temp = merge_quantification_annotation(df_quantification, df_annotation,
                                              tab_level=tab_level+1)
    df_temp[sample_name] = df_temp[sample_name].astype(int)
    df_temp["start"] = df_temp["start"].astype(int)
    df_temp["end"] = df_temp["end"].astype(int)


    # 获取单个样本的gene信息
    gene_info = load_annotation(input_filename=annotation_location,
                                additional_info=annotation_col,
                                only_type="gene",
                                tab_level=tab_level+1)
    expressed_gene = list(set(df_temp["annot_gene_id"].to_list()))
    expressed_gene = pandas.DataFrame(expressed_gene)
    gene_info = pandas.merge(expressed_gene, gene_info, left_on=0, right_on="gene_id", how="left")
    gene_info["start"] = gene_info["start"].astype(int)
    gene_info["end"] = gene_info["end"].astype(int)

    # 记录单个样本的gene信息
    exist_gene_id = {}  # 存储重复出现的gene的id {<new_gene_id>: existed_gene_id}
    for i in gene_info.index:
        chr = gene_info.at[i, "chr"]
        gene_id = gene_info.at[i, "gene_id"]
        gene_name = gene_info.at[i, "gene_name"]
        source = gene_info.at[i, "source"]
        strand = gene_info.at[i, "strand"]
        start = gene_info.at[i, "start"]
        end = gene_info.at[i, "end"]

        # gene_mark有三种值: "True", "False", <gene_id> 此处的<gene_id>代表与当前Gene对象重复的gene的id
        gene_mark = total.add_gene(Gene(chr=chr, gene_id=gene_id, gene_name=gene_name, gene_biotype="un_classified",
                                        source=source, strand=strand, start=start, end=end,
                                        transcript_dict={}, exon_dict={}))
        if gene_mark == "True" or gene_mark == "False":
            pass
        else:
            exist_gene_id[gene_id] = gene_mark


    # 获取单个样本的exon信息
    exon_info = load_annotation(input_filename=annotation_location,
                                additional_info=annotation_col,
                                only_type="exon",
                                tab_level=tab_level+1)
    expressed_gene = list(set(df_temp["annot_transcript_id"].to_list()))
    expressed_gene = pandas.DataFrame(expressed_gene)
    exon_info = pandas.merge(expressed_gene, exon_info, left_on=0, right_on="transcript_id")
    exon_info["start"] = exon_info["start"].astype(int)
    exon_info["end"] = exon_info["end"].astype(int)

    # 记录单个样本的exon信息
    exon_info_group = exon_info.groupby("transcript_id")
    for transcript_id in exon_info_group.groups.keys():
        temp_df = exon_info_group.get_group(transcript_id)

        gene_id = temp_df.at[temp_df.index[0], "gene_id"]
        # 由于单个样本的gene信息已在先前步骤中被统计，所以，此时再对同一样本的gene添加exon信息时,
            # 不会出现未被记录的gene. 也就是说，此处exon所对应的gene可分为两种情况：
            #   1.该gene_id已被记录  2.该gene_id未被记录但已记录相应的gene信息(即gene信息相同但id不同)
        if gene_id in exist_gene_id.keys():
            # 该exon所在的gene属于重复出现的gene(该gene_id被映射到其他的gene_id)
            gene_id = exist_gene_id.get(gene_id)
        else:
            # 该exon所在的gene未被映射到其他gene_id
            pass

        # 向指定gene_id添加exon信息
        for i in temp_df.index:
            total.gene_dict[gene_id].add_exon(temp_df.at[i,"start"], temp_df.at[i,"end"])


    # 记录单个样本的transcript信息
    # 根据df_temp
    # 向每一个gene添加相关的transcript_id, transcript的start及end, 在单个样本中的counts数量
    temp_transcript_info = df_temp.set_index("annot_transcript_id")
    exon_info_group = exon_info.groupby("transcript_id")
    for transcript_id in temp_transcript_info.index:
        temp_exon_info = exon_info_group.get_group(transcript_id)

        transcript_name = temp_transcript_info.at[transcript_id, "annot_transcript_name"]
        start = temp_transcript_info.at[transcript_id, "start"]
        end = temp_transcript_info.at[transcript_id, "end"]
        sample_name = sample_name
        sample_counts = temp_transcript_info.at[transcript_id, sample_name]
        exon_list = [[temp_exon_info.at[i, "start"], temp_exon_info.at[i, "end"]] for i in temp_exon_info.index]

        gene_id = temp_transcript_info.at[transcript_id, "annot_gene_id"]
        if gene_id in exist_gene_id.keys():
            gene_id = exist_gene_id.get(gene_id)
        else:
            pass

        total.gene_dict[gene_id].add_transcript(transcript_id=transcript_id,
                                                transcript_name=transcript_name,
                                                transcript_biotype="un_classified",
                                                start=start, end=end,
                                                exon_list=exon_list,
                                                sample_name=sample_name,
                                                sample_counts=sample_counts)


    # 获取单个样本的sample_name并添加到total的对象中
    sample_name = sample_name
    sample_list = total.sample_list
    sample_list.append(sample_name)
    total.sample_list = sample_list

    return total

@log
def save_df(total, file_path, tab_level=0):
    """
    input:
        total, Total, 存储数据的对象
        file_path, str, 要保存的文件的绝对路径
    change:
        将Total对象转为pandas.DataFrame并保存到指定位置
    output:
        df, pandas.DataFrame, 由Total转换而来的df
    """
    total = total
    file_path = file_path

    df = total.get_df(total.sample_list)
    df = df.sort_values(by="gene_id")
    df = df.sort_values(by="chr")
    df.to_csv(file_path, sep='\t', index=None)

    return df


@log
def get_biotype_info(annotation_file_location, tab_level=0):
    """
    input:
        annotation_file_location, str, the path of annotation file
    change:
        get biotype of gene or transcript from annotation file
    output:
        biotype_dict, dict, {<gene_id>: <gene_biotype>, <gene_id>: <gene_biotype>, ...
                             <transcript_id>: <transcript_biotype>, <transcript_id>: <transcript_biotype>, ...}
    """
    annotation_file_location = annotation_file_location

    # load annotation file
    if annotation_file_location[-3:] == ".gz":
        df = gzip.open(annotation_file_location, 'rb')
    else:
        df = open(annotation_file_location, 'r')

    # get biotype info of gene/transcript
    biotype_dict = {}
    for line in df:
        temp_dict = {}

        line = line.decode()
        line = line.split('\t')

        # only retain row of gene/transcript
        if len(line) <= 2:
            continue
        else:
            if line[2] == "gene":
                gene_mark = True
                transcript_mark = False
            elif line[2] == "transcript":
                gene_mark = False
                transcript_mark = True
            else:
                continue

        # get info of biotype
        line = line[-1]
        line = line.strip('\n')
        line = [i.strip(' ') for i in line.split(';')]

        if '' in line:
            line.remove('')

        for i in line:
            i = i.replace('"', '')
            i = i.split(' ')
            temp_dict[i[0]] = i[1]

        if gene_mark is True:
            biotype_dict[temp_dict.get("gene_id")] = temp_dict.get("gene_biotype")
        elif transcript_mark is True:
            biotype_dict[temp_dict.get("transcript_id")] = temp_dict.get("transcript_biotype")
        else:
            df.close()
            raise ValueError("[Error]type of line is {} of:\n\t {}\n".format(line[2], line))
    df.close()

    return biotype_dict


def add_biotype_info(total, biotype_dict, tab_level=0):
    """
    input:
        total, Total, the object of Total
        biotype_dict, dict, {<gene_id>: <gene_biotype>, ... , <transcript_id>: <transcript_biotype>, ...}
    change:
        1.遍历total中的每一个gene及相应的transcripts
        2.去掉gene_id/transcript_id的版本号
        3.寻找biotype_dict中是否有对应的biotype信息, 若有则向gene/transcript添加相应的biotype, 若无则更改值为"un_classified"
        4.打印具有相应biotype的gene数及transcript数'
    output:
        bool, return True if perform successed else return False
    """
    total = total
    biotype_dict = biotype_dict

    num_gene = 0
    num_transcript = 0

    for gene_id, gene_object in total.gene_dict.items():
        # 只修改已有注释信息的gene或transcript
        if gene_id[0:4] != "ENSG":
            pass
        else:
            num_gene +=1
            temp_gene_id = gene_id.split('.')[0]
            total.gene_dict[gene_id].gene_biotype = biotype_dict.get(temp_gene_id, "un_classified")

        # 检验该gene所包含的每一个transcript_id开头是否为"ENST", 若是则修改transcript_biotype, 若否则跳过该transcript
        for transcript_id in gene_object.transcript_dict.keys():
            if transcript_id[0:4] != "ENST":
                continue
            else:
                num_transcript +=1
                temp_transcript_id = transcript_id.split('.')[0]
                total.gene_dict[gene_id].transcript_dict[transcript_id]["transcript_biotype"] = biotype_dict.get(temp_transcript_id, "un_classified")
    print("\t[result]the counts of gene updated biotype: {}".format(num_gene))
    print("\t[result]the counts of transcript updated biotype: {}".format(num_transcript))

    return True


@log
def get_cell_line_sample_dict(sample_info, GEO_index=False, tab_level=0):
    """
    input:
        sample_info, df, the df of sample_info.tsv
        GEO_index, bool, GEO_id is index
    change:
        以dict形式存储cellLine与样本之间的对应关系
    output:
        cell_line_sample_dict, dict, {disease1: [sample1, sample2,...], ...}
    """
    sample_info = sample_info
    GEO_index = GEO_index

    if GEO_index is False:
        pass
    else:
        sample_info = sample_info.reset_index()

    cell_line_sample_dict = {}
    sample_info_group = sample_info.groupby("cell_line")
    for cellLine in sample_info_group.groups.keys():
        cell_line_sample_dict[cellLine] = sample_info_group.get_group(cellLine)["GEO_accession"].to_list()

    print("{}[Result]cell_line_sample_dict: {}".format('\t'*(tab_level+1), cell_line_sample_dict))

    return cell_line_sample_dict


@log
def load_sample_info(filename, tab_level=0):
    """
    input:
        读取sample_info.tsv文件
    change:
        读取sample_info文件, 对其中的disease列, 以下划线替换空格
    output:
        sample_info, pandas.DataFrame
    """
    filename = filename

    sample_info = pandas.read_csv(filename, sep='\t')
    sample_info["disease"] = sample_info["disease"].map(lambda x: x.replace(' ', '_'))

    return sample_info


def dedupExonList(exonList):
    """
    input:
        exonList, list, [[<exonStart>, <exonEnd>], [<exonStart>, <exonEnd>], ...]
    change:
        在exonList中, 各exon的范围存在重叠
        该函数遍历整个exonList, 合并存在重叠的exon
    output:
        dedupExonList, list, 去重叠后的exon的list
    """
    if True:
        exonList = exonList

    # 改变exonList中的exon, 使exonStart <= exonEnd
    exonList = list(map(lambda x: [x[1], x[0]] if x[0]>x[1] else [x[0], x[1]], exonList))

    # 按照exonStart的大小进行升序排序
    exonList = sorted(exonList, key=lambda x: x[0])

    dedupExonList = []
    [exonStart, exonEnd] = exonList[0]
    for [tempExonStart, tempExonEnd] in exonList:
        if exonStart <= tempExonStart <= exonEnd:
            # 当前exon与上一个exon存在重叠
            exonEnd = max(exonEnd, tempExonEnd)
        else:
            # 当前exon与上一个exon不重叠
            dedupExonList.append([exonStart, exonEnd])
            exonStart = tempExonStart
            exonEnd = tempExonEnd
    dedupExonList.append([exonStart, exonEnd])

    return dedupExonList
