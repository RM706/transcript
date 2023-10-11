import os
import re
import pandas
import time

__functionVersion__ = "V1.0(Editor) 2023-10-08"

'''
函数说明:
    log()                           装饰器
    读取ENCODE数据用函数
        get_file_path()             一级函数    记录ENCODE的annotation及quantification文件绝对路径
        读取quantificaton文件
            loadQuantification()    一级函数    读取ENCODE的quantification数据
        读取annotation文件
            exonIndexMapped()       一级函数    修改exonIndex中可被映射的exonId
            pickEncodeAnnotation()  一级函数    格式化处理ENCODE基因组注释中的每一行信息
            loadAnnotation()        二级函数    读取单个annotation文件中exon|transcript|gene的数据
        loadSample()                三级函数    读取单个样本的quantificaton及annotation数据
        loadData()                  四级函数    读取所有样本的数据
    读取参考基因组注释用函数
        pickRefAnnotation()         一级函数    读取参考基因组注释文件数据
'''

# 修饰器
def log(function):
    name = function.__name__

    def wrapper(*args, **kwargs):
        tab = kwargs.get("tabLevel", 0)

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


# 一级函数
def get_file_path(dataPath, sampleNameList, tabLevel=0):
    """
    input:
        dataPath, str, data文件夹路径
        sample_name_list, list, 含有样本名的列表
        tabLevel, int, the parament for log format, there is no need to change
    change:
        搜索file_path下的每个样本文件夹中的文件名,
            记录annotation及quantification文件的绝对路径
    output:
        sample_file_location, dict, 该字典存储每个样本的annotation及quantification文件的绝对路径
                                    {<sample1>:{"annotation": <dir1>, "quantification": <dir2>},
                                     <sample2>:{"annotation": <dir1>, "quantification": <dir2>}, ...}
    """

    dataPath = dataPath
    sampleNameList = sampleNameList

    sample_file_location = {sample: {} for sample in sampleNameList}
    for sample in sampleNameList:
        filename_list = os.listdir("{}{}".format(dataPath, sample))
        for filename in filename_list:
            if re.search(pattern="annotation", string=filename) is not None:
                sample_file_location[sample]["annotation"] = "{}{}/{}".format(dataPath, sample, filename)
            elif re.search(pattern="quantification", string=filename) is not None:
                sample_file_location[sample]["quantification"] = "{}{}/{}".format(dataPath, sample, filename)
            else:
                continue

    return sample_file_location


# 一级函数
def pickEncodeAnnotation(info, allowType=None, tabLevel=0):
    '''
    input:
        info,\t str,\t annotation中的一行信息\n
        allowType,\t list,\t =["gene", "transcript", "exon"], 需要进行分析的行的类型\n
    change:
        格式化处理基因组注释中的每一行信息\n
    output:
        info,\t dict,\t 以键值对形式存储了info中的信息\n
                        若info不属于分析的范围, 则返回None\n
    '''
    info = info
    allowType = (allowType, ["gene", "transcript", "exon"])[allowType is None]

    # 过滤掉注释的行
    if '#' in info:
        return None

    info = info.split('\t')

    chr = info[0]
    lineType = info[2]
    start = info[3]
    end = info[4]
    strand = info[6]

    # 不处理allowType中不存在的类型
    if lineType not in allowType:
        return None

    # 处理annotation文件中最后一列的信息
    info = info[-1]
    info = info.replace('\n', '')
    info = info.split(';')
    (None, info.remove(''))['' in info]  # 去除列表中的''
    info = [i.strip(' ') for i in info]  # 去除最前面的空格
    info = [i.split(' ') for i in info]  # 分离key与value
    # 转换为dict格式
    info = {i[0]: i[1:] for i in info}
    for k, v in info.items():
        v = ''.join(v)
        v = v.replace('"', '')
        info[k] = v

    # 向dict中添加gene的基础信息
    info["chr"] = chr
    info["type"] = lineType
    info["start"] = start
    info["end"] = end
    info["strand"] = strand

    return info


# 一级函数
def pickRefAnnotation(info, allowType=None, tabLevel=0):
    '''
    input:
        info,\t str,\t annotation中的一行信息\n
        allowType,\t list,\t =["gene", "transcript", "exon"], 需要进行分析的行的类型\n
    change:
        格式化处理参考基因组注释中的每一行信息\n
    output:
        info,\t dict,\t 以键值对形式存储了info中的信息\n
                        若info不属于分析的范围, 则返回None\n
    '''
    info = info
    allowType = (allowType, ["gene", "transcript", "exon"])[allowType is None]

    # 检查参数格式
    if not isinstance(allowType, list):
        raise ValueError("pickRefAnnotation() the allowType show be list")

    # 过滤掉注释的行
    if '#' in info:
        return None

    info = info.split('\t')

    chr = info[0]
    lineType = info[2]
    start = info[3]
    end = info[4]
    strand = info[6]

    # 不处理allowType中不存在的类型
    if lineType not in allowType:
        return None

    # 处理annotation文件中最后一列的信息
    info = info[-1]
    info = info.replace('\n', '')
    info = info.split(';')
    (None, info.remove(''))['' in info]  # 去除列表中的''
    info = [i.strip(' ') for i in info]  # 去除最前面的空格
    info = [i.split(' ') for i in info]  # 分离key与value
    # 转换为dict格式
    info = {i[0]: i[1:] for i in info}
    for k, v in info.items():
        v = ''.join(v)
        v = v.replace('"', '')
        info[k] = v

    # 向dict中添加gene的基础信息
    info["chr"] = chr
    info["type"] = lineType
    info["start"] = start
    info["end"] = end
    info["strand"] = strand

    # 将版本号添加到id中
    if "exon" in allowType:
        info["exon_id"] = "{}.{}".format(info["exon_id"], info["exon_version"])
    if "transcript" in allowType:
        info["transcript_id"] = "{}.{}".format(info["transcript_id"], info["transcript_version"])
    if "gene" in allowType:
        info["gene_id"] = "{}.{}".format(info["gene_id"], info["gene_version"])

    return info


# 一级函数
def loadQuantification(filename, cutoff):
    '''
    input:
        filename, str
        cutoff, int
    change:
        读取quantification文件, 并根据CUTOFF进行过滤, 筛选得到 rep值>=cutoff值 的transcript
    '''
    filename = filename
    cutoff = cutoff

    df = pandas.read_csv(filename, sep='\t')
    df = df.loc[df[df.columns[-1]]>=cutoff, ["annot_gene_id", "annot_transcript_id", df.columns[-1]]]
    df = df.rename(columns={df.columns[-1]: "rep", "annot_gene_id": "geneId", "annot_transcript_id": "transcriptId"})

    return df


# 一级函数
def exonIndexMapped(total, exonIndex):
    '''
    input:
        total, Total Object
        exonIndex, dict, {<chr>: {<transcriptId>: [exon1, exon2, exon3, ...], ...}, ...}
    change:
        根据total.exonExisted修改exonIndex中的exonId, 将可以映射的exonId映射到existedExonId
    output:
        exonIndex, dict
    '''
    total = total
    exonIndex = exonIndex

    for chr, transcriptIdDict in exonIndex.items():
        for transcriptId, exonList in transcriptIdDict.items():
            temp = exonList
            marker = False  # 是否需要修改exonList
            for i in range(0, len(exonList)):
                exonId = temp[i]
                # 如果exonId可映射到其他exonId
                if exonId in total.exonExisted.keys():
                    marker = True
                    existedExonId = total.exonExisted[exonId]
                    temp[i] = existedExonId
            if marker is True:
                exonIndex[chr][transcriptId] = temp

    return exonIndex


# 二级函数
def loadAnnotation(total, sample, filename, type, geneIdSet, transcriptIdSet, chrList, exonIndex, countsTranscript):
    '''
    input:
        total, Total
        sample, str
        filename, str, ENCODE的annotation文件的位置
        type, str, "gene" or "transcript" or "exon"
        geneIdSet, set, 只有存在于该set中的gene才会被纳入分析
        transcriptIdSet, set, 只有存在于该set中的transcript才会被纳入分析
        chrList, list, 存储有标准染色体 ["chrX", "chrY", "chrM", "chr1", "chr2", ...]
        exonIndex, dict, 存储有transcript对应的exon信息{<chr>: {<transcriptId>: [<exon1>, <exon2>, <exon3>, ...]}, ...}
        countsTranscript, dict, 存储有transcript在sample中的counts数{<transcriptId>: <counts>}
    change:
        读取ENCODE的annotation文件中指定类型的数据至变量total中
    output:
        list, [total, None或者dfExon]已加载annotation数据的total
    '''
    total = total
    sample = sample
    filename = filename
    type = type
    geneIdSet = geneIdSet
    transcriptIdSet = transcriptIdSet
    chrList = chrList
    exonIndex = exonIndex
    countsTranscript = countsTranscript

    # 检查参数是否符合格式
    if type not in ["gene", "transcript", "exon"]:
        raise ValueError("[Error]loadAnnotation() the parament type is {}, not is 'gene' or 'transcript' or 'exon'".format(type))

    with open(filename, 'r') as file:
        if type == "gene":
            for line in file:
                line = pickEncodeAnnotation(info=line, allowType=[type])
                if line != None:
                    geneId = line.get("gene_id", "")
                    geneName = line.get("gene_name", "")
                    geneStatus = line.get("gene_status", "")
                    geneBiotype = line.get("gene_type", "")
                    chr = line.get("chr", "")
                    start = int(line.get("start", ""))
                    end = int(line.get("end", ""))
                    strand = line.get("strand", "")
                    # 过滤掉位于非标准染色体上的gene
                    if chr not in chrList:
                        continue
                    else:
                        pass
                    if geneId in geneIdSet:
                        total.geneAdd(status=geneStatus, chr=chr, strand=strand, start=start, end=end, geneId=geneId, geneName=geneName, geneBiotype=geneBiotype)
                    else:
                        continue
            return [total, None]
        elif type == "transcript":
            for line in file:
                line = pickEncodeAnnotation(info=line, allowType=[type])
                if line is None:
                    continue
                geneId = line.get("gene_id", "")
                transcriptId = line.get("transcript_id", "")
                transcriptName = line.get("transcript_name", "")
                transcriptStatus = line.get("transcript_status", "")
                transcriptBiotype = line.get("transcript_type", "")
                chr = line.get("chr", "")
                start = int(line.get("start", ""))
                end = int(line.get("end", ""))
                # 过滤掉位于非标准染色体上的transcript
                if chr not in chrList:
                    continue
                else:
                    pass
                # 过滤掉不存在于transcriptIdSet中的exonId, 以去除噪音
                if transcriptId not in transcriptIdSet:
                    continue
                # 检查该geneId是否被映射到了其他的geneId
                geneId = total.geneIdGet(geneId=geneId)
                # 添加transcript信息
                total.geneDict[geneId].transcriptAdd(status=transcriptStatus,
                                                     transcriptId=transcriptId,
                                                     transcriptName=transcriptName,
                                                     transcriptBiotype=transcriptBiotype,
                                                     start=start,
                                                     end=end,
                                                     exonList=exonIndex[chr][transcriptId],
                                                     sample=sample,
                                                     counts=countsTranscript[transcriptId])
            return [total, None]
        else:
            exonIndex = {i: {} for i in chrList}  # {<chr>: {<transcriptId>: [<exon1>, <exon2>, ...], ...}, ...}
            for line in file:
                line = pickEncodeAnnotation(info=line, allowType=["exon"])
                if line != None:
                    geneId = line.get("gene_id", "")
                    transcriptId = line.get("transcript_id", "")
                    exonId = line.get("exon_id", "")
                    exonStatus = line.get("exon_status", "")
                    chr = line.get("chr", "")
                    strand = line.get("strand", "")
                    start = int(line.get("start", ""))
                    end = int(line.get("end", ""))
                    # 过滤掉位于非标准染色体上的exon
                    if chr not in chrList:
                        continue
                    else:
                        pass
                    # 过滤掉不存在于transcriptIdSet中的exonId, 以去除噪音
                    if transcriptId not in transcriptIdSet:
                        continue
                    # 检查exonId是否为纯数字
                    if "ENSE" not in exonId:
                        # 如果exonId为纯数字, 则在exonId前添加样本名
                        exonId = "{}_{}".format(sample, exonId)
                    total.exonAdd(geneId=geneId, exonId=exonId, status=exonStatus, chr=chr, strand=strand, start=start, end=end)
                    # 此时的indexExon中存在有可被映射的exonId, 需要在后续步骤中进行处理
                    temp = exonIndex[chr].get(transcriptId, [])
                    temp.append(exonId)
                    exonIndex[chr][transcriptId] = temp

            return [total, exonIndex]


# 三级函数
def loadSample(total, sampleCellline, sampleFilePath, sample, cutoff, chrList):
    '''
    input:
        total, Total
        sampleCellline, dict, 含有sample所对应的cellline信息{<sample>: <cellline>, ...}
        sampleFilePath, dict, 含有每个sample对应的annotation及quantification文件路径
        sample, str, 样本名
        cutoff, int, quantification文件中只有rep值>=该值的transcript才会被纳入分析
        chrList, list, 存储有标准染色体 ["chrX", "chrY", "chrM", "chr1", "chr2", ...]
    change:
        加载对应样本的annotation及quantification数据
    '''
    total = total
    sampleCellline = sampleCellline
    sampleFilePath = sampleFilePath
    sample = sample
    cutoff = cutoff
    chrList = chrList

    fileQuantification = sampleFilePath[sample]["quantification"]
    fileAnnotation = sampleFilePath[sample]["annotation"]

    # 加载quantification数据
    dfQuantification = loadQuantification(filename=fileQuantification, cutoff=cutoff)
    geneIdSet = set(dfQuantification["geneId"])
    transcriptIdSet = set(dfQuantification["transcriptId"])
    ## 获取transcript的counts数
    dfQuantification = dfQuantification.loc[:, ["transcriptId", "rep"]]
    dfQuantification = dfQuantification.set_index("transcriptId")
    countsTranscript = dfQuantification.to_dict()["rep"]

    # 加载annotation数据
    ## 加载gene数据
    total = loadAnnotation(total=total, sample=sample, filename=fileAnnotation, type="gene", geneIdSet=geneIdSet, transcriptIdSet=transcriptIdSet, chrList=chrList, exonIndex=None, countsTranscript=None)[0]
    ## 加载exon数据
    [total, exonIndex] = loadAnnotation(total=total, sample=sample, filename=fileAnnotation, type="exon", geneIdSet=geneIdSet, transcriptIdSet=transcriptIdSet, chrList=chrList, exonIndex=None, countsTranscript=None)
    exonIndex = exonIndexMapped(total=total, exonIndex=exonIndex)
    ## 加载transcript数据
    ### 过滤掉位于非标准染色体上的transcript
    transcriptIdSet = set()
    for exonDict in exonIndex.values():
        transcriptIdSet = transcriptIdSet.union(set(exonDict.keys()))
    countsTranscript = {k: v for k, v in countsTranscript.items() if k in transcriptIdSet}
    ### 添加transcript信息
    total = loadAnnotation(total=total, sample=sample, filename=fileAnnotation, type="transcript", geneIdSet=geneIdSet, transcriptIdSet=transcriptIdSet, chrList=chrList, exonIndex=exonIndex, countsTranscript=countsTranscript)[0]

    # 补充sample信息
    cellline = sampleCellline[sample]
    temp = total.celllineInfo.get(cellline, [])
    temp.append(sample)
    total.celllineInfo[cellline] = temp

    return total


# 四级函数
def loadData(total, sampleCellline, sampleFilePath, cutoff, chrList, tabLevel=0):
    '''
    change:
        读取所有样本的数据
    '''
    total = total
    sampleCellline = sampleCellline
    sampleFilePath = sampleFilePath
    cutoff = cutoff
    chrList = chrList

    progress = 0
    progressTotal = len(sampleFilePath.keys())
    for sample in sampleFilePath.keys():
        progress += 1
        print("{}[Progress]{}/{}".format('\t'*(tabLevel+1), progress, progressTotal), end='\r')
        total = loadSample(total=total, sampleCellline=sampleCellline, sampleFilePath=sampleFilePath, sample=sample, cutoff=cutoff, chrList=chrList)

    return total
