import numpy

classVersion = "V1.5(Editor) 2023-10-17"

'''
提示:
    在对exon, transcript, gene的对象进行增, 删, 合并后, 需要进行的操作
        refresh()
        reIndex()
        computeRelativeExpression()
        computeCellLineExpression()
'''

'''
函数说明:
    杂项
        strSplit()          一级函数    将含有version的id分割为[<id>, <version>]
'''

# add思路
'''
0. 分别声明total, gene, transcript, exon的数据格式
    0.0 Total
        0.0.0 属性
            0.0.0.0 geneDict, dict, 存储Gene对象
            0.0.0.1 exonDict, dict, 存储Exon对象
            0.0.0.2 celllineInfo, dict, {<cellLine1>: [<sample1>, <sample2>, ...], ...}
            0.0.0.3 geneIndex, dict, 根据位置建立的对gene的索引{{<chr>: '+': {<start>: {<end>: <gene id>,...},...},'-': {...},...},}
            0.0.0.4 geneExisted, dict, 记录重复的gene id之间的映射关系
        0.0.1 方法
            0.0.1.0 __geneCheck, 检查gene是否已被记录
            0.0.1.1 geneAdd, 向geneDict中添加Gene对象
            0.0.1.2 dfGeneGet, 将Total中的所有gene信息整理成df并返回
            0.0.1.3 dfTranscriptGet, 将Total中的所有transcript信息整理成df并返回
            0.0.1.4 dfExonGet, 将Total中的所有exon信息整理成df并返回
            0.0.1.5 computeRelativeExpression, 计算transcript在每个样本中的相对表达量
            0.0.1.6 computeCellLineExpression, 计算transcript在每个细胞系中的相对表达量
    0.1 Gene
        0.1.0 属性
            0.1.0.0 status, str, 
            0.1.0.1 chr, str, gene所在的染色体号
            0.1.0.2 strand, str, gene的链的方向
            0.1.0.3 start, int, gene的start的位置(与strand无关, start<end)
            0.1.0.4 end, int, gene的end的位置(与strand无关, start<end)
            0.1.0.5 geneId, str, 带有版本号
            0.1.0.6 geneName, str
            0.1.0.7 geneBiotype, str
            0.1.0.9 exonList, list, 存储该gene所包含的所有exon
            0.1.0.10 transcriptDict, dict, 存储映射的Transcript对象, key为id
        0.1.1 方法
            0.1.1.0 geneGet, 以dict格式获取gene的信息
            0.1.1.1 __exonCheck, 检查exon是否已被记录
            0.1.1.2 exonAdd, 向exonList中添加Exon对象
            0.1.1.3 __transcriptCheck, 检查transcript是否已被记录
            0.1.1.4 transcriptAdd, 向transcriptDict中添加Transcript对象
    0.2 Exon
        0.2.0 属性
            0.2.0.0 status, str
            0.2.0.1 exonId, str
            0.2.0.2 exonStart, int, (与strand无关, start<end)
            0.2.0.3 exonEnd, int, (与strand无关, start<end)
        0.2.1 方法
            0.2.1.0 exonGet, 以dict格式获取exon的信息
    0.3 Transcript
        0.3.0 属性
            0.3.0.0 status, str
            0.3.0.1 transcriptId, str
            0.3.0.2 transcriptName, str
            0.3.0.3 transcriptBiotype, str
            0.3.0.4 transcriptStart, int, (与strand无关, start<end)
            0.3.0.5 transcriptEnd, int, (与strand无关, start<end)
            0.3.0.6 exonList, list, 该转录本的exon组成
            0.3.0.7 countsExpression, dict, {<sample_name1>: <expression>, <sample_name2>: <expression>, ...}
            0.3.0.8 relativeExpression, dict, {<sample_name1>: <expression>, <sample_name2>: <expression>, ...}
            0.3.0.9 cellLineExpression, dict, {<cellLine1>: <expression>, <cellLine2>: <expression>, ...}
        0.3.1 方法
            0.3.1.0 transcriptGet, 以dict格式获取transcript的信息

1. 逐个样本读取quantification文件及annotation文件
    1.0 读取quantification文件数据
        1.0.0 设定阈值COUNTSCUTOFF, 只有rep值>=该阈值的transcript才会被纳入接下来的分析
        1.0.1 声明set类型变量geneIdSet, 存储过滤后的transcript
        1.0.2 声明set类型变量transcriptIdSet, 存储过滤后的transcript
    1.1 读取annotation文件数据
        1.1.1 添加gene信息
            声明变量
                声明一个dict类型变量geneExisted, 记录重复的gene id之间的映射关系
                声明一个dict类型变量geneIndex, 记录gene的位置信息  {{<chr>: '+': {<start>: {<end>: <gene id>,
                                                                                                    ...},
                                                                                            ...},
                                                                                    '-': {...},
                                                                            ...},
                                                                            }
            逐行读取annotation中的gene数据
            过滤掉不存在于geneIdSet中的geneId, 以去除噪音
            检查gene id是否已记录
                检查gene id是否已映射到其他gene id
                    若已映射, 则不处理该gene
                    若未映射, 则进行下一步
                检查gene id是否已记录
                    若gene id已记录, 则不分析该gene
                    若gene id未记录, 则检查是否存在这样一种情况: gene已被记录但对应到了其他的gene id
                        检查该gene的start是否已存在于geneIndex[<chr>][<strand>].keys()中
                            若不存在, 则认为该gene未被记录, 向total中添加该gene对象, 并向geneIndex添加位置信息
                            若存在, 则检查该gene的end是否已存在于geneIndex[<chr>][<strand>][<start>].keys()
                                若不存在, 则认为该gene未被记录, 向total中添加该gene对象, 并向geneIndex添加位置信息
                                若存在, 则认为该gene已被记录, 向geneExisted中添加gene id的映射关系 {<new gene id>: <existed gene id>}
        1.1.2 添加exon信息
            逐行读取annotation中的exon数据
            声明变量
                声明一个dict类型变量exonExisted, 记录重复的exon id之间的映射关系{<new exon id>: <existed exon id>, ...}
                声明一个dict类型变量exonIndex, 记录exon的位置关系{<chr>: {<start>: {<end>: <exon id>,
                                                                                            ...},
                                                                                    ...},
                                                                            ...}
            过滤掉不存在于transcriptIdSet中的exonId, 以去除噪音
            检查exonId是否为纯数字
                如果exonId为纯数字, 则在exonId前添加样本名
            检查该exon是否位于标准染色体上
                若不位于标准染色体上, 则不分析该exon
                若位于标准染色体上, 则进行下一步
            检查该exon是否已被记录
                检查该exonId是否已映射到其他exonId
                    若已映射, 则检查对应的gene中是否已记录<existed exon id>
                        若已记录, 则不处理该exon
                        若未记录, 则在对应的gene中添加该exonId
                    若未映射, 则进行下一步
                检查exonId是否已存在于total.exonDict.keys()
                    若存在, 则检查对应gene中是否已记录该exon
                        若存在, 则忽略该exon
                        若不存在, 则在对应gene中记录该exon
                    若不存在, 则检查该exon的start是否存在于exonIndex[<chr>].keys()
                        若不存在, 则认为该exon未被记录
                            向相应的gene添加该exonId
                            向exonDict中添加该exon信息
                            向exonIndex中添加该exon信息
                        若存在, 则检查该exon的end是否存在于exonIndex[<chr>][<start>].keys()中
                            若不存在, 则认为该exon未被记录
                                向相应的gene添加该exonId
                                向exonDict中添加该exon信息
                                向exonIndex中添加该exon信息
                            若存在, 则认为该exon已被记录只是exonId不同
                                向exonExisted中添加该映射{<old exon id>: <existed exon id>}
                                检查对应gene中是否已记录<existed exon id>
                                    若未记录, 则向相应的gene添加<existed exon id>
        1.1.3 添加transcript信息
            声明一个dict类型变量transcriptExisted, 记录重复的transcript id之间的映射关系
            逐行读取annotation中的transcript数据
            检查该transcriptId是否映射到了其他transcriptId
                若已映射, 则认为该transcript已存在
                    将该transcript的{<sample>: <counts>}添加到total.geneDict[<geneId>].transcriptDict[<existed transcriptId>].countsExpressoin中
                    结束处理
                若未映射, 则进行下一步
            检查在transcript对应的gene中该transcriptId是否已存在
                若存在, 则认为该transcript已记录
                    在total.geneDict[<geneId>].transcriptDict[<transcript id>].countsExpression中添加{<sample>: <counts>}
                    结束处理
                若不存在, 则检查transcript的start是否存在于transcriptIndex.keys()
                    若不存在, 则transcript未记录
                        在total.geneDict[<geneId>].transcriptDict中添加transcript对象
                        在total.geneDict[<geneId>].transcriptDict[transcriptId].countsExpression中添加{<sample>: <counts>}
                        在total.geneDict[<geneId>].transcriptIndex中添加transcript位置信息
                    若存在, 则检查transcript的end是否存在于transcriptIndex[<start>].keys()
                        若不存在, 则transcript未记录
                            在total.geneDict[<geneId>].transcriptDict中添加transcript对象
                            在total.geneDict[<geneId>].transcriptDict[transcriptId].countsExpression中添加{<sample>: <counts>}
                            在total.geneDict[<geneId>].transcriptIndex中添加transcript位置信息
                        若存在, 则检查transcript的exon组成是否与相应的transcript相同
                            遍历映射的transcript, 获取每个transcript的exon组成
                            若存在transcript的exon组成与该transcript的exon组成一致,则认为该transcript已被记录
                                在total.geneDict[<geneId>].transcriptExisted中添加映射信息{<new transcript id>: <existed transcript id>}
                                在total.geneDict[<geneId>].transcriptDict[<existed transcript id>].countsExpression中添加{<sample>: <counts>}
                            若不存在, 则认为该transcript未被记录
                                在total.geneDict[<geneId>].transcriptDict中添加transcript对象
                                在total.geneDict[<geneId>].transcriptDict[transcriptId].countsExpression中添加{<sample>: <counts>}
                                在total.geneDict[<geneId>].transcriptIndex中添加transcript位置信息
'''


# 一级函数
def strSplit(str):
    '''
    change:
        将str按照<split>进行分割, 形成一个list
            如果无法分割, 就新建一个list=[<str>, -1]
    '''
    if '.' in str:
        str = str.split('.')
    else:
        str = [str, -1]
    return str


# 定义Total
class Total(object):
    '''
    属性
        geneDict, dict, 存储Gene对象
        geneIndex, dict, 根据位置建立的对gene的索引{<chr>: {<strand>: {<start>: {<end>: <gene id>, ...}, ...}, ...}, ...}
        geneExisted, dict, 记录重复的gene id之间的映射关系
        exonDict, dict, 存储Exon对象
        exonIndex, dict, 根据位置建立的对exon的索引
        exonExisted, dict, 记录重复的exon id之间的映射关系
        celllineInfo, dict, {<cellline1>: [<sample1>, <sample2>, ...], ...}
    方法
        _exonCheck()                检查exon是否已记录
        _exonExistedAdd()           向exonExisted中添加映射
        _exonExistedRefresh()       去除exonExisted中的多重映射
        _exonMerge()                合并两个exon
        _transcriptMerge()          合并两个transcript, 可跨gene合并
        _geneCheck()                检查gene是否已记录
        _geneExistedAdd()           向geneExisted中添加映射
        _geneMerge()                合并两个gene
        _addTSSTES()                自检并更新gene, transcript, exon的TSS及TES信息
        exonAdd()                   添加exon
        exonUseGet()                获取使用指定exon的所有gene及transcript
        exonIndexBuild()            建立符合过滤条件的exon的index
        transcriptIdGet()           输入一个transcriptId或transcriptName, 返回一个已考虑映射关系的transcriptId或None
        transcriptIndexBuild()      建立符合过滤条件的transcript的index
        geneAdd()                   添加gene
        geneIdGet()                 输入一个geneId或geneName, 返回一个已考虑映射关系的geneId或None
        geneIndexBuild()            建立符合过滤条件的gene的index
        refresh()                   自检所有gene及transcript的exonList中是否存在可映射的exon, 若存在, 则进行映射
        reIndex()                   重新建立self.geneIndex, self.exonIndex及每个geneObject中的transcriptIndex
        computeRelativeExpression() 计算每个transcript在每个样本中的相对表达量, 并统计gene及exon在每个样本中的相对表达量
        computeCellLineExpression() 计算每个gene, transcript, exon在每个细胞系中的相对表达量
    '''
    def __init__(self, geneDict=None, geneIndex=None, geneExisted=None, exonDict=None, exonIndex=None, exonExisted=None, celllineInfo=None, chrList=None):
        if chrList is None:
            chrList = ["chr{}".format(i) for i in range(1,22+1)]
            chrList.append("chrX")
            chrList.append("chrY")
            chrList.append("chrM")
        self.geneDict = (geneDict, {})[geneDict is None]
        self.geneIndex = (geneIndex, {i: {'+': {}, '-': {}} for i in chrList})[geneIndex is None]
        self.geneExisted = (geneExisted, {})[geneExisted is None]
        self.exonDict = (exonDict, {})[exonDict is None]
        self.exonIndex = (exonIndex, {i: {} for i in chrList})[exonIndex is None]
        self.exonExisted = (exonExisted, {})[exonExisted is None]
        self.celllineInfo = (celllineInfo, {})[celllineInfo is None]

    def _exonCheck(self, exonId, chr, start, end):
        '''
        change:
            检查该exon是否已被记录
        output:
            list, 若exonId已记录, 则返回[0, None], 旧exonId
                  若exonId未记录 且start未记录, 则返回[1, None], 新exonId
                  若exonId未记录 且start已记录 且end未记录, 则返回[2, None], 新exonId
                  若exonId未记录 且start已记录 且end已记录, 则返回[3, <existed exon id>], 旧exonId
        '''
        exonId = exonId
        chr = chr
        start = start
        end = end

        # 检查该exon是否已被记录
        if exonId in self.exonDict.keys():
            # exonId已记录
            return [0, None]
        else:
            # exonId未记录,
            # 检查该exon的start是否已被记录
            temp = self.exonIndex[chr]
            if start not in temp.keys():
                # exonId未记录 且start未记录
                return [1, None]
            else:
                # exonId未记录 且start已记录
                # 检查该exon的end是否已被记录
                temp = temp[start]
                if end not in temp.keys():
                    # exonId未记录 且start已记录 且end未记录
                    return [2, None]
                else:
                    # exonId未记录 且start已记录 且end已记录
                    existedExonId = temp[end]
                    return [3, existedExonId]

    def _exonExistedAdd(self, oldExonId, newExonId):
        '''
        change:
            向exonExisted中添加映射
        return:
            bool, 正常情况下返回True, 但若这个映射是自我映射则返回False
        '''
        oldExonId = oldExonId
        newExonId = newExonId

        keys = set(self.exonExisted.keys())
        values = set(self.exonExisted.values())

        # 寻找到所有可以映射到newExonId的exonId
        oldExonIdSet = {oldExonId}
        workSet = {oldExonId}
        while workSet:
            exonId = workSet.pop()
            if exonId in values:
                # 有其他的exonId映射到了该exonId
                temp = {k for k, v in self.exonExisted.items() if v==exonId}
                # workSet添加这些exonId, 检查是否还有exonId可映射到这些exonId        
                workSet = workSet.union(temp)
            else:
                oldExonIdSet.add(exonId)
        # 寻找到最终的newExonId
        while newExonId in keys:
            newExonId = self.exonExisted[newExonId]

        # 不进行自我映射
        if oldExonId == newExonId:
            return False
        elif newExonId in oldExonIdSet:
            return False

        for exonId in oldExonIdSet:
            self.exonExisted[exonId] = newExonId

        return True

    def _exonExistedRefresh(self):
        '''
        change:
            一般情况下应该用不上这个函数
            修改exonExisted中的A-->B-->C-->D这种映射
                A-->D
                B-->D
                C-->D
        '''
        temp = set()

        keys = set(self.exonExisted.keys())
        values = set(self.exonExisted.values())
        for k, v in self.exonExisted.items():
            if k in values or v in keys:
                temp.add((k, v))

        while temp:
            (k, v) = temp.pop()
            self._exonExistedAdd(oldExonId=k, newExonId=v)

        return None

    def _exonMerge(self, exonId1, exonId2):
        '''
        change:
            注意: exon在merge之后可能会导致transcript产生重合, 还应在执行该函数后检查transcript是否存在重复
            合并两个exon
                不可靠exon --> 可靠exon
                倾向于保留exonId1
            0. 分别检查exon1及exon2是否已记录
            1. 确当映射方向并映射
                若两个exon的status为"KNOWN", 则在self.exonExisted中将exon2映射到exon1
                若一个exon的status为"KNOWN", 则在self.exonExisted中将非KNOWN的exon映射到KNOWN的exon
                若没有exon的status为"KNOWN", 则在self.exonExisted中将exon2映射到exon1
            2. 在self.exonDict中删除不可靠exon
            3. 返回一个dict, 保存有{deletedExonId: exonId}
        return:
            tuple, 正常情况下返回(<deletedExonId>: <exonId>), 但当deletedExonId==remainedExonId时会返回("", "")
        '''
        exonId1 = exonId1
        exonId2 = exonId2

        deletedExonId = None
        remainedExonId = None
        # 分别检查exon1及exon2是否已记录
        if exonId1 not in self.exonDict.keys():
            raise ValueError("[Error]exonMerge() the exonId1 {} not in self.exonDict".format(exonId1))
        elif exonId2 not in self.exonDict.keys():
            raise ValueError("[Error]exonMerge() the exonId2 {} not in self.exonDict".format(exonId2))
        elif exonId1 == exonId2:
            print("[Warning]_exonMerge()--exonId1 {} == exonId2 {}".format(exonId1, exonId2))
            return {"": ""}
        else:
            pass

        # 检查exon1及exon2的status
        status1 = self.exonDict[exonId1].status
        status2 = self.exonDict[exonId2].status
        if status1 == "KNOWN" and status2 == "KNOWN":
            # 将exon2映射到exon1
            deletedExonId = exonId2
            remainedExonId = exonId1
        elif status1 == "KNOWN":
            # 将exon2映射到exon1
            deletedExonId = exonId2
            remainedExonId = exonId1
        elif status2 == "KNOWN":
            # 将exon1映射到exon2
            deletedExonId = exonId1
            remainedExonId = exonId2
        else:
            # 将exon2映射到exon1
            deletedExonId = exonId2
            remainedExonId = exonId1
        marker = self._exonExistedAdd(oldExonId=deletedExonId, newExonId=remainedExonId)

        # 处理异常: deletedExonId == remainedExonId
        if marker is False:
            return ("", "")

        # 修改self.exonIndex
        try:
            self.exonIndex[self.exonDict[deletedExonId].chr][self.exonDict[deletedExonId].start][self.exonDict[deletedExonId].end] = remainedExonId
        except KeyError as k:
            # 目前出现问题的可能原因是: 在exon经过校正后其start及end发生了改变
            #   但在self.exonIndex中其start及end的位置并没有发生改变
            #   所以在根据exon的新的start/end修改self.exonIndex时, self.exonIndex中无法找到对应的位置
            pass
            #print(k)
            #print("[Warning]_exonMerge()--remained:{}, delete:{}, start:{}, end:{}".format(remainedExonId, deletedExonId, self.exonDict[deletedExonId].start, self.exonDict[deletedExonId].end))

        # self.exonDict中删除不可靠exon
        self.exonDict.pop(deletedExonId)

        return (deletedExonId, remainedExonId)

    def _transcriptMerge(self, transcriptId1, transcriptId2, geneId1, geneId2):
        '''
        input:
            transcriptId1, str
            transcriptId2, str
            geneId1, str, this value should be geneId of transcriptId1
            geneId2, str, this value should be geneId of transcriptId2
        change:
            合并两个transcript, 倾向于保留transcriptId1
                0. 检查两transcriptId是否位于同一gene对象中
                    若位于同一gene对象中, 则直接调用gene对象中的_transcriptMerge()
                1. 若两transcriptId位于不同的gene
                    检查两transcriptId是否可映射, 若可映射则进行映射
                    检查两transcript是否存在, 若存在则进行下一步
                    确定映射方向并映射
                        若2个transcript的status为KNOWN, 则在self.transcriptExisted中将transcriptId2映射到transcriptId1
                        若1个transcript的status为KNOWN, 则在self.transcriptExisted中将非KNOWN的transcriptId映射到KNOWN的transcriptId
                        若0个transcript的status为KNOWN, 则在self.transcriptExisted中将transcriptId2映射到transcriptId1
                    更改remainedTranscript的countsExpression, exonList
                    删除transcriptDict中deletedTranscript对象
                    更新两gene的exonList  refresh()
                    重建两gene的transcriptIndex  reIndex()
        return:
            tuple, 正常情况下返回(<deletedTranscriptId>, <remainedTranscriptId>), 当deletedTranscriptId==remainedTrarnscriptId时会返回("", "")
        '''
        transcriptId1 = transcriptId1
        transcriptId2 = transcriptId2
        geneId1 = geneId1
        geneId2 = geneId2

        if geneId1 == geneId2:
            # 两transcript位于同一gene
            # 调用Gene对象中的_transcriptMerge()
            (deletedTranscriptId, remainedTranscriptId) = self.geneDict[geneId1]._transcriptMerge(transcriptId1=transcriptId1, transcriptId2=transcriptId2)
            return (deletedTranscriptId, remainedTranscriptId)
        else:
            # 两transcript位于不同gene
            # 检查两transcriptId是否可映射
            if transcriptId1 in self.geneDict[geneId1].transcriptExisted.keys():
                transcriptId1 = self.geneDict[geneId1].transcriptExisted[transcriptId1]
            if transcriptId2 in self.geneDict[geneId2].transcriptExisted.keys():
                transcriptId2 = self.geneDict[geneId2].transcriptExisted[transcriptId2]
            # 检查两transcript是否存在, 若存在则进行下一步
            if transcriptId1 not in self.geneDict[geneId1].transcriptDict.keys():
                raise ValueError("_transcriptMerge(): transcriptId1 {} not in geneId1 {}".format(transcriptId1, geneId1))
            if transcriptId2 not in self.geneDict[geneId2].transcriptDict.keys():
                raise ValueError("_transcriptMerge(): transcriptId2 {} not in geneId2 {}".format(transcriptId2, geneId2))
            # 确定映射方向并映射
            transcriptObject1 = self.geneDict[geneId1].transcriptDict[transcriptId1]
            transcriptObject2 = self.geneDict[geneId2].transcriptDict[transcriptId2]
            if transcriptObject1.status=="KNOWN" and transcriptObject2.status=="KNOWN":
                # 将transcriptId2映射到transcriptId1
                deletedTranscriptId = transcriptId2
                deletedGeneId = geneId2
                remainedTranscriptId = transcriptId1
                remainedGeneId = geneId1
            elif transcriptObject1.status=="KNOWN":
                # 将transcriptId2映射到transcriptId1
                deletedTranscriptId = transcriptId2
                deletedGeneId = geneId2
                remainedTranscriptId = transcriptId1
                remainedGeneId = geneId1
            elif transcriptObject2.status=="KNOWN":
                # 将transcriptId1映射到transcriptId2
                deletedTranscriptId = transcriptId1
                deletedGeneId = geneId1
                remainedTranscriptId = transcriptId2
                remainedGeneId = geneId2
            else:
                # 将transcriptId2映射到transcriptId1
                deletedTranscriptId = transcriptId2
                deletedGeneId = geneId2
                remainedTranscriptId = transcriptId1
                remainedGeneId = geneId1
            marker = self.geneDict[remainedGeneId]._transcriptExistedAdd(oldTranscriptId=deletedTranscriptId, newTranscriptId=remainedTranscriptId)
            # 处理异常: deletedTranscriptId == remainedTranscriptId
            if marker is False:
                return ("", "")
            # 更改remainedTranscript的countsExpression
            for sample, counts in self.geneDict[deletedGeneId].transcriptDict[deletedTranscriptId].countsExpression.items():
                temp = self.geneDict[remainedGeneId].transcriptDict[remainedTranscriptId].countsExpression.get(sample, 0)
                temp = temp + counts
                self.geneDict[remainedGeneId].transcriptDict[remainedTranscriptId].countsExpression[sample] = temp
            # 更改remainedTranscript的exonList
            deletedExon = self.geneDict[deletedGeneId].transcriptDict[deletedTranscriptId].exonList.copy()
            remainedExon = self.geneDict[remainedGeneId].transcriptDict[remainedTranscriptId].exonList.copy()
            remainedExon = remainedExon + deletedExon
            remainedExon = list(set(remainedExon))
            self.geneDict[remainedGeneId].transcriptDict[remainedTranscriptId].exonList = remainedExon
            # 删除transcriptDict中deletedTranscript对象
            self.geneDict[deletedGeneId].transcriptDict.pop(deletedTranscriptId)
            # 更新两gene的exonList
            self.geneDict[deletedGeneId].refresh(exonExisted=self.exonExisted)
            self.geneDict[remainedGeneId].refresh(exonExisted=self.exonExisted)
            # 重建两gene的transcriptIndex
            self.geneDict[deletedGeneId].reIndex()
            self.geneDict[remainedGeneId].reIndex()

            return (deletedTranscriptId, remainedTranscriptId)

    def _geneCheck(self, geneId, chr, strand, start, end):
        '''
        input:
            geneId, str
            chr, str
            strand, str
            start, int
            end, int
        change:
            检查gene是否已被记录
        output:
            list, 当 geneId已记录 时, 返回[0, None], 旧gene
                当 geneID未记录 且start未记录 时, 返回[1, None], 新gene
                当 geneID未记录 且start已记录 且end未记录 时, 返回[2, None], 新gene
                当 geneId未记录 且start已记录 且end已记录 时, 返回[3, <existed gene id>], 旧gene
        '''
        # 检查gene id是否已记录
        if geneId in self.geneDict.keys():
            # gene id已记录, 旧gene
            return [0, None]
        else:
            # gene id未记录
            # 检查该gene的start是否已记录
            temp = self.geneIndex[chr][strand]
            if start not in temp.keys():
                # 该gene的start未记录, 新gene
                return [1, None]
            else:
                # 该gene的start已记录
                # 检查该gene的end是否已记录
                temp = temp[start]
                if end not in temp.keys():
                    # 该gene的end未被记录, 新gene
                    return [2, None]
                else:
                    # 该gene的end已被记录, 旧gene
                    return [3, temp[end]]

    def _geneExistedAdd(self, oldGeneId, newGeneId):
        '''
        change:
            向geneExisted中添加映射
        return:
            bool, True--已进行映射, False--该映射为自我映射
        '''
        oldGeneId = oldGeneId
        newGeneId = newGeneId

        keys = set(self.geneExisted.keys())
        values = set(self.geneExisted.values())

        # 寻找到所有可以映射到newGeneId的geneId
        oldGeneIdSet = {oldGeneId}
        workSet = {oldGeneId}
        while workSet:
            geneId = workSet.pop()
            if geneId in values:
                # 有其他的geneId映射到了该geneId
                temp = {k for k, v in self.geneExisted.items() if v==geneId}
                # workSet添加这些geneId, 检查是否还有geneId可映射到这些geneId        
                workSet = workSet.union(temp)
            else:
                oldGeneIdSet.add(geneId)
        # 寻找到最终的newGeneId
        while newGeneId in keys:
            newGeneId = self.geneExisted[newGeneId]

        # 不进行自我映射
        if oldGeneId == newGeneId:
            return False
        elif newGeneId in oldGeneIdSet:
            return False

        for geneId in oldGeneIdSet:
            self.geneExisted[geneId] = newGeneId

        return True

    def _geneMerge(self, geneId1, geneId2):
        '''
        change:
            合并两gene
            0. 检验两gene是否存在
            1. 确定映射方向并映射
                若2个gene的status为KNOWN, 则在self.geneExisted中将geneId2映射到geneId1
                若1个gene的status为KNOWN, 则在self.geneExisted中将非KNOWN的geneId映射到KNOWN的geneId
                若0个gene的status为KNOWN, 则在self.geneExisted中将geneId2映射到geneId1
            2. 合并gene
                合并exonList
                合并transcriptExisted
                    直接update
                合并gene中的每个transcript
                    transcriptId不同可直接添加, 在重建index过程中就可以合并重复的transcript
                    transcriptId相同则需要合并countsExpression
                重建transcriptIndex
            3. 删除deletedGene
        return:
            tuple, 正常情况下返回(<deletedGeneId>, <remainedGeneId>), 但当deletedExonId==remainedExonId时会返回("", "")
        '''
        geneId1 = geneId1
        geneId2 = geneId2

        if geneId1 not in self.geneDict.keys():
            raise ValueError("[Error]_geneMerge() the geneId1 {} not in self.geneDict".format(geneId1))
        elif geneId2 not in self.geneDict.keys():
            raise ValueError("[Error]_geneMerge() the geneId2 {} not in self.geneDict".format(geneId2))
        else:
            pass

        # 确定映射方向
        geneObject1 = self.geneDict[geneId1]
        geneObject2 = self.geneDict[geneId2]
        if geneObject1.status=="KNOWN" and geneObject2.status=="KNOWN":
            # 将geneId2映射到geneId1
            deletedGeneId = geneId2
            remainedGeneId = geneId1
        elif geneObject1.status=="KNOWN":
            # 将geneId2映射到geneId1
            deletedGeneId = geneId2
            remainedGeneId = geneId1
        elif geneObject2.status=="KNOWN":
            # 将geneId1映射到geneId2
            deletedGeneId = geneId1
            remainedGeneId = geneId2
        else:
            # 将geneId2映射到geneId1
            deletedGeneId = geneId2
            remainedGeneId = geneId1
        marker = self._geneExistedAdd(oldGeneId=deletedGeneId, newGeneId=remainedGeneId)

        # 处理异常: deletedGeneId == remainedGeneId
        if marker is False:
            return ("", "")

        # 合并exonList
        temp = self.geneDict[remainedGeneId].exonList + self.geneDict[deletedGeneId].exonList
        temp = list(set(temp))  # 去重
        self.geneDict[remainedGeneId].exonList = temp

        # 合并transcriptExisted
        self.geneDict[remainedGeneId].transcriptExisted.update(self.geneDict[deletedGeneId].transcriptExisted)

        # 合并gene中的每个transcript
        for transcriptId, transcriptObject in self.geneDict[deletedGeneId].transcriptDict.items():
            # 判断transcript的id是否已记录
            if transcriptId in self.geneDict[remainedGeneId].transcriptDict.keys():
                # transcriptId已记录于remainedGene中
                self._transcriptMerge(transcriptId1=transcriptId, geneId1=remainedGeneId,
                                      transcriptId2=transcriptId, geneId2=deletedGeneId)
            else:
                # transcriptId未记录于remainedGene中
                # 可直接添加
                # 注意, 此时可能会在gene中引入重复的transcript, 但是在后续的重建index过程中可以去除
                temp = {transcriptId: transcriptObject}
                self.geneDict[remainedGeneId].transcriptDict.update(temp)

        # 重建transcriptIndex
        self.geneDict[remainedGeneId].reIndex()

        # 修改geneIndex
        try:
            self.geneIndex[self.geneDict[deletedGeneId].chr][self.geneDict[deletedGeneId].start][self.geneDict[deletedGeneId].end] = remainedGeneId
        except KeyError as k:
            # 目前出现问题的可能原因是: 在gene经过校正后其start及end发生了改变
            #   但在self.geneIndex中其start及end的位置并没有发生改变
            #   所以在根据gene的新的start/end修改self.geneIndex时, self.geneIndex中无法找到对应的位置
            pass

        # 删除deletedGene
        self.geneDict.pop(deletedGeneId)

        return (deletedGeneId, remainedGeneId)

    def _addTSSTES(self):
        for geneObject in self.geneDict.values():
            geneObject._addTSSTES()

        for exonObject in self.exonDict.values():
            exonObject._addTSSTES()

        return None

    def exonAdd(self, geneId, exonId, exonVersion, status, chr, strand, start, end):
        '''
        change:
            添加exon信息
        '''
        geneId = geneId
        exonId = exonId
        exonVersion = exonVersion
        status = status
        chr = chr
        strand = strand
        start = start
        end = end

        # 检查exon是否位于标准染色体列表中
        if chr not in self.exonIndex.keys():
            return None

        # 检查对应的geneId是否被映射到了其他的geneId
        geneId = self.geneIdGet(geneId=geneId)

        # 检查该exon是否已被记录
        ## 检查该exonId是否已映射到其他exonId
        if exonId in self.exonExisted.keys():
            existedExonId = self.exonExisted[exonId]
            # 若已映射, 则检查对应的gene中是否已记录<existed exon id>
            if self.geneDict[geneId]._exonCheck(existedExonId) is True:
                # 若已记录, 则不处理该exon
                return None
            else:
                # 若未记录, 则在对应的gene中添加该exonId
                self.geneDict[geneId].exonList.append(existedExonId)
                return None
        else:
            pass
        ## 检查该exon是否可映射到其他exonId
        checkMarker = self._exonCheck(exonId=exonId, chr=chr, start=start, end=end)
        if checkMarker[0] == 0:
            # exonId已记录
            marker = self.geneDict[geneId]._exonCheck(exonId)
            if marker is True:
                # 对应的gene已记录该exonId
                return None
            else:
                # 向对应的gene的exonList中添加该exonId
                self.geneDict[geneId].exonList.append(exonId)
                return None
        elif checkMarker[0] == 1:
            # exonId未记录 且start未记录
            self.geneDict[geneId].exonList.append(exonId)
            self.exonDict[exonId] = Exon(status=status, exonId=exonId, exonVersion=exonVersion, chr=chr, strand=strand, start=start, end=end)
            self.exonIndex[chr][start] = {end: exonId}
        elif checkMarker[0] == 2:
            # exonId未记录 且start已记录 且end未记录
            self.geneDict[geneId].exonList.append(exonId)
            self.exonDict[exonId] = Exon(status=status, exonId=exonId, exonVersion=exonVersion, chr=chr, strand=strand, start=start, end=end)
            self.exonIndex[chr][start][end] = exonId
        elif checkMarker[0] == 3:
            # exonId未记录 且start已记录 且end已记录
            existedExonId = checkMarker[1]
            # 判断exon的映射方向
            existedExonStatus = self.exonDict[existedExonId].status
            if status == existedExonStatus == "KNOWN":
                # 两exon都为KNOWN, 新exon-->旧exon
                self._exonExistedAdd(oldExonId=exonId, newExonId=existedExonId)
                if existedExonId not in self.geneDict[geneId].exonList:
                    self.geneDict[geneId].exonList.append(existedExonId)
            elif status == "KNOWN":
                # 新的exon为KNOWN, 旧exon-->新exon
                self.exonDict[exonId] = Exon(status=status, exonId=exonId, exonVersion=exonVersion, chr=chr, strand=strand, start=start, end=end)
                self._exonMerge(exonId1=exonId, exonId2=existedExonId)
                if exonId not in self.geneDict[geneId].exonList:
                    self.geneDict[geneId].exonList.append(exonId)
            elif existedExonStatus == "KNOWN":
                # 旧的exon为KNOWN, 新exon-->旧exon
                self._exonExistedAdd(oldExonId=exonId, newExonId=existedExonId)
                if existedExonId not in self.geneDict[geneId].exonList:
                    self.geneDict[geneId].exonList.append(existedExonId)
            else:
                # 两exon都不为KNOWN, 新exon-->旧exon
                self._exonExistedAdd(oldExonId=exonId, newExonId=existedExonId)
                if existedExonId not in self.geneDict[geneId].exonList:
                    self.geneDict[geneId].exonList.append(existedExonId)
        else:
            raise ValueError("exonAdd() -- the chechMarker is {}, only 0,1,2,3".format(checkMarker))

    def exonUseGet(self, exonId):
        '''
        change:
            遍历所有gene, 寻找使用指定exonId的gene
        return:
            {<geneId>: {transcriptId, ...}, ...}
        '''
        exonId = exonId

        result = {}
        for geneId, geneObject in self.geneDict.items():
            exonInGene = geneObject.exonUseGet(exonId=exonId)
            result.update(exonInGene)

        return result

    def exonIndexBuild(self, siteType="start", statusSet=None):
        '''
        input:
            siteType, str, "start" or "end" or "TSS" or "TES"
            statusSet, set, 只有status在这个set中的exon才会被纳入index
        return:
            index, dict, {<chr>: {<siteType>: {<siteType2>: <exonId>}, ...}, ...}
        '''
        siteType = siteType
        statusSet = (statusSet, None)[statusSet is None]

        mapped = {"start": "end", "end": "start", "TSS": "TES", "TES": "TSS"}
        index = {}
        for exonId, exonObject in self.exonDict.items():
            chr = exonObject.chr
            status = exonObject.status
            if statusSet is not None:
                # 过滤掉不符合条件的gene
                if status not in statusSet:
                    continue
            data = {"start": exonObject.start, "end": exonObject.end, "TSS": exonObject.TSS, "TES": exonObject.TES}
            if chr not in index.keys():
                index[chr] = {data[siteType]: {data[mapped[siteType]]: exonId}}
            else:
                if data[siteType] not in index[chr].keys():
                    index[chr][data[siteType]] = {data[mapped[siteType]]: exonId}
                else:
                    if data[mapped[siteType]] not in index[chr][data[siteType]].keys():
                        index[chr][data[siteType]][data[mapped[siteType]]] = exonId
                    else:
                        pass

        return index

    def transcriptIdGet(self, transcriptId=None, transcriptName=None):
        '''
        change:
            遍历所有gene, 寻找是否有对应的transcriptId或transcriptName
        return:
            dict or bool, 若能找到transcriptId, 则返回{"geneId": <geneId>, "transcriptId": <transcriptId>}
                          若未能找到transcriptId, 则返回None
        '''
        transcriptId = transcriptId
        transcriptName = transcriptName

        if transcriptId is None and transcriptName is None:
            raise ValueError("transcriptIdGet() need transcriptId or transcriptName")
        elif transcriptName is None:
            # 根据transcriptId进行寻找
            for geneId, geneObject in self.geneDict.items():
                result = geneObject.transcriptIdGet(transcriptId=transcriptId)
                if result is not None:
                    # transcriptId已找到
                    return {"geneId": geneId, "transcriptId": result}
            # transcriptId未找到
            return None
        else:
            # 根据transcriptName进行寻找
            for geneId, geneObject in self.geneDict.items():
                result = geneObject.transcriptIdGet(transcriptName=transcriptName)
                if result is not None:
                    # transcriptName已找到
                    return {"geneId": geneId, "transcriptId": result}
            # transcriptName未找到
            return None

    def transcriptIndexBuild(self, siteType="start", statusSet=None):
        '''
        input:
            siteType, str, "start" or "end" or "TSS" or "TES"
            statusSet, set, 只有status在这个set中的transcript才会被纳入index
        return:
            index, dict, {<chr>: {<siteType>: {siteType2: geneId}, ...}, ...}
        '''
        siteType = siteType
        statusSet = (statusSet, None)[statusSet is None]

        mapped = {"start": "end", "end": "start", "TSS": "TES", "TES": "TSS"}
        index = {}

        for geneId, geneObject in self.geneDict.items():
            chr = geneObject.chr
            for transcriptId, transcriptObject in geneObject.transcriptDict.items():
                status = transcriptObject.status
                # 过滤掉不符合条件的gene
                if statusSet is not None:
                    if status not in statusSet:
                        continue
                data = {"start": transcriptObject.start, "end": transcriptObject.end, "TSS": transcriptObject.TSS, "TES": transcriptObject.TES}
                if chr not in index.keys():
                    index[chr] = {data[siteType]: {data[mapped[siteType]]: [transcriptId]}}
                else:
                    if data[siteType] not in index[chr].keys():
                        index[chr][data[siteType]] = {data[mapped[siteType]]: [transcriptId]}
                    else:
                        if data[mapped[siteType]] not in index[chr][data[siteType]].keys():
                            index[chr][data[siteType]][data[mapped[siteType]]] = [transcriptId]
                        else:
                            index[chr][data[siteType]][data[mapped[siteType]]].append(transcriptId)

        return index

    def geneAdd(self, status, chr, strand, start, end, geneId, geneVersion, geneName, geneBiotype):
        '''
        input:
            status, str
            chr, str
            strand, str
            start, int
            end, int
            geneId, str
            geneVersion, int
            geneName, str
            geneBiotype, str
        change:
            添加Gene对象
        '''
        status = status
        chr = chr
        strand = strand
        start = start
        end = end
        geneId = geneId
        geneVersion = geneVersion
        geneName = geneName
        geneBiotype = geneBiotype

        # 检查gene是否位于标准染色体列表中
        if chr not in self.geneIndex.keys():
            return None

        # 检查gene是否已记录
        ## 检查gene id是否已映射到其他gene id
        if geneId in self.geneExisted.keys():
            # gene id已映射到其他gene id
            return None
        else:
            # gene id未映射到其他gene id
            pass
        ## 检查gene id是否已记录
        checkMarker = self._geneCheck(geneId=geneId, chr=chr, strand=strand, start=start, end=end)
        if checkMarker[0] == 0:
            # geneId已记录, 旧gene
            return None
        elif checkMarker[0] == 1:
            # geneID未记录 且start未记录, 新gene
            self.geneDict[geneId] = Gene(status=status,
                                         chr=chr,
                                         strand=strand,
                                         start=start,
                                         end=end,
                                         geneId=geneId,
                                         geneVersion=geneVersion,
                                         geneName=geneName,
                                         geneBiotype=geneBiotype)
            self.geneIndex[chr][strand][start] = {end: geneId}
            return None
        elif checkMarker[0] == 2:
            # geneID未记录 且start已记录 且end未记录, 新gene
            self.geneDict[geneId] = Gene(status=status,
                                         chr=chr,
                                         strand=strand,
                                         start=start,
                                         end=end,
                                         geneId=geneId,
                                         geneVersion=geneVersion,
                                         geneName=geneName,
                                         geneBiotype=geneBiotype)
            self.geneIndex[chr][strand][start][end] = geneId
            return None
        elif checkMarker[0] == 3:
            # geneId未记录 且start已记录 且end已记录, 旧gene
            existedGeneId = checkMarker[1]
            # TODO:确定gene的映射方向 当前会出现KNOWN的gene映射到NOVEL的gene
            self._geneExistedAdd(oldGeneId=geneId, newGeneId=existedGeneId)
            return None
        else:
            raise ValueError("[Error]geneAdd() -- checkMarker is {}".format(checkMarker))

    def geneIdGet(self, geneId=None, geneName=None):
        '''
        input:
            geneId, str, =None
            geneName, str, =None
        change:
            如果输入geneId, 则获取geneDict中的geneId, 已考虑映射情况
            如果数据geneName, 则根据输入的geneName返回该geneName所对应的geneId
        output:
            str or bool, 如果寻找得到geneId则返回geneId, 否则返回None
        '''
        geneId = geneId
        geneName = geneName

        if geneId is None and geneName is None:
            raise ValueError("geneIdGet() need geneId or geneName")
        elif geneName is None:
            # 根据geneId进行寻找
            if geneId in self.geneExisted.keys():
                # geneId被映射到了其他geneId
                geneId = self.geneExisted[geneId]
                return geneId
            else:
                # geneId未映射到其他geneId
                if geneId in self.geneDict.keys():
                    # 该geneId已被记录
                    return geneId
                else:
                    # 该geneId未被记录
                    return None
        else:
            # 根据geneName进行寻找
            for geneId, geneObject in self.geneDict.items():
                if geneName == geneObject.geneName:
                    # 已根据geneName寻找到对应的geneId
                    return geneId
            # 未能根据geneName寻找到geneId
            return None

    def geneIndexBuild(self, siteType="start", statusSet=None):
        '''
        input:
            siteType, str, "start" or "end" or "TSS" or "TES"
            statusSet, set, 只有status在这个set中的gene才会被纳入index
        return:
            index, dict, {<chr>: {<siteType>: {siteType2: geneId}, ...}, ...}
        '''
        siteType = siteType
        statusSet = (statusSet, None)[statusSet is None]

        mapped = {"start": "end", "end": "start", "TSS": "TES", "TES": "TSS"}
        index = {}
        for geneId, geneObject in self.geneDict.items():
            chr = geneObject.chr
            status = geneObject.status
            # 过滤掉不符合条件的gene
            if statusSet is not None:
                if status not in statusSet:
                    continue
            data = {"start": geneObject.start, "end": geneObject.end, "TSS": geneObject.TSS, "TES": geneObject.TES}
            if chr not in index.keys():
                index[chr] = {data[siteType]: {data[mapped[siteType]]: geneId}}
            else:
                if data[siteType] not in index[chr].keys():
                    index[chr][data[siteType]] = {data[mapped[siteType]]: geneId}
                else:
                    if data[mapped[siteType]] not in index[chr][data[siteType]].keys():
                        index[chr][data[siteType]][data[mapped[siteType]]] = geneId
                    else:
                        pass

        return index

    def refresh(self):
        '''
        change:
            对所有gene及所有transcript的exonList进行自检, 将可映射的exon进行映射(当exon的映射发生改变后需要进行该操作)
        return:
            dict, {<geneId>: {"geneChanged": <bool>, "transcriptList": [<transcriptId1>, <transcriptId2>, ...]}, ...}
        '''

        changedDict = {}  # 存储exonCombination发生了改变的gene及transcript
        for geneId, geneObject in self.geneDict.items():
            [geneChangedMarker, changedTranscriptIdList] = geneObject.refresh(exonExisted=self.exonExisted)
            if geneChangedMarker is False and len(changedTranscriptIdList)==0:
                # gene及该gene下的所有transcript的exonCombination均未发生变化
                pass
            elif geneChangedMarker is False:
                # gene的exonCombination未改变, 但存在transcript的exonCombination发生改变
                changedDict[geneId] = {"geneChanged": False, "transcriptList": changedTranscriptIdList}
            elif len(changedTranscriptIdList) == 0:
                # gene的exonCombination发生改变, transcript的exonCombination均未发生改变
                changedDict[geneId] = {"geneChanged": True, "transcriptList": changedTranscriptIdList}
            else:
                # gene的exonCombination发生改变, transcript的exonCombination也发生了改变
                changedDict[geneId] = {"geneChanged": True, "transcriptList": changedTranscriptIdList}

        return changedDict

    def reIndex(self):
        '''
        change:
            重新建立self.geneIndex, self.exonIndex及每个geneObject中的transcriptIndex
            合并gene
                合并exonList
                合并transcriptExisted
                    直接update
                合并transcriptDict
                    transcriptId不同可直接添加, 在重建index过程中就可以合并重复的transcript
                    transcriptId相同则需要合并countsExpression
                重建transcriptIndex
        '''

        # 重新建立self.exonIndex
        newIndex = {}
        tempExonDict = self.exonDict.copy()
        for exonId, exonObject in tempExonDict.items():
            chr = exonObject.chr
            start = exonObject.start
            end = exonObject.end

            if chr not in newIndex.keys():
                # 该exon的chr未记录
                newIndex[chr] = {start: {end: exonId}}
            else:
                # 该exon的chr已记录
                if start not in newIndex[chr].keys():
                    # 该exon的start未记录
                    newIndex[chr][start] = {end: exonId}
                else:
                    # 该exon的start已记录
                    if end not in newIndex[chr][start].keys():
                        # 该exon的end未记录
                        newIndex[chr][start][end] = exonId
                    else:
                        # 该exon的end已记录, 表明该exon已重复
                        existedExonId = newIndex[chr][start][end]
                        (deletedExonId, remainedExonId) = self._exonMerge(exonId1=existedExonId, exonId2=exonId)
                        newIndex[chr][start][end] = remainedExonId
        self.exonIndex = newIndex

        # 重新建立transcriptIndex
        ## 先将transcript中可以进行映射的exon进行映射, 因为经过上一步后exonExisted可能会发生改变
        for geneId, geneObject in self.geneDict.items():
            geneObject.refresh(exonExisted=self.exonExisted)
        ## 再重新在每个gene中建立transcriptIndex
        for geneId, geneObject in self.geneDict.items():
            geneObject.reIndex()

        # 重新建立self.geneIndex
        newIndex = {}
        tempGeneDict = self.geneDict.copy()
        for geneId, geneObject in tempGeneDict.items():
            chr = geneObject.chr
            strand = geneObject.strand
            start = geneObject.start
            end = geneObject.end
            if chr not in newIndex.keys():
                # 该gene的chr未记录
                newIndex[chr] = {strand: {start: {end: geneId}}}
            else:
                # 该gene的chr已记录
                if strand not in newIndex[chr].keys():
                    # 该gene的strand未记录
                    newIndex[chr][strand] = {start: {end: geneId}}
                else:
                    # 该gene的strand已记录
                    if start not in newIndex[chr][strand].keys():
                        # 该gene的start未记录
                        newIndex[chr][strand][start] = {end: geneId}
                    else:
                        # 该gene的start已记录
                        if end not in newIndex[chr][strand][start].keys():
                            # 该gene的end未记录
                            newIndex[chr][strand][start][end] = geneId
                        else:
                            # 该gene的end已记录
                            # 认为两gene可以合并, 尽量保留existedGene
                            existedGeneId = newIndex[chr][strand][start][end]
                            (deletedGeneId, remainedGeneId) = self._geneMerge(geneId1=existedGeneId, geneId2=geneId)
                            newIndex[chr][strand][start][end] = remainedGeneId
        self.geneIndex = newIndex

        # 重新建立transcriptIndex
        for geneId, geneObject in self.geneDict.items():
            geneObject.reIndex()

        # 更新gene, transcript, exon的TSS及TES
        self._addTSSTES()

        return None

    def computeRelativeExpression(self):
        '''
        change:
            计算每个transcript在每个样本中的相对表达量, 并统计gene及exon在每个样本中的相对表达量
                0. 遍历self.celllineInfo中的values, 获取所有样本名
                1. 遍历所有transcript, 统计该样本中所有样本的总的counts数
                2. 处理transcript
                    对每一个样本, 计算transcript在该样本中的相对表达量
                        相对表达量 = 在该样本中该transcript的counts * 10^6 / 该样本中所有transcript的counts数的总和
                3. 处理gene
                    根据gene所包含的transcript, 统计该gene的counts数
                    计算gene在每个样本中的相对表达量
                4. 处理exon
                    根据transcript所包含的exon, 统计exon的counts数
                    计算exon在每个样本中的相对表达量
        '''

        # 获取所有样本名
        sampleDict = {}  # {<sample>: <total counts>, ...}
        for sampleList in self.celllineInfo.values():
            for sample in sampleList:
                sampleDict[sample] = None

        # 统计所有样本的总的counts数
        for sample in sampleDict.keys():
            totalCounts = 0
            for geneId, geneObject in self.geneDict.items():
                for transcriptId, transcriptObject in geneObject.transcriptDict.items():
                    totalCounts = totalCounts + transcriptObject.countsExpression.get(sample, 0)
            sampleDict[sample] = totalCounts

        # 处理transcript
        # 计算每个transcript的相对表达量
        for sample in sampleDict.keys():
            for geneId, geneObject in self.geneDict.items():
                for transcriptId, transcriptObject in geneObject.transcriptDict.items():
                    relativeExpression = transcriptObject.countsExpression.get(sample, 0) * 10**6 / sampleDict[sample]
                    transcriptObject.relativeExpression[sample] = relativeExpression

        # 处理gene
        # 统计每个gene的counts
        for geneObject in self.geneDict.values():
            for sample in sampleDict.keys():
                counts = 0
                for transcriptObject in geneObject.transcriptDict.values():
                    counts = counts + transcriptObject.countsExpression.get(sample, 0)
                geneObject.countsExpression[sample] = counts
        # 计算每个gene的relativeExpression
        for geneObject in self.geneDict.values():
            for sample, counts in geneObject.countsExpression.items():
                relativeExpression = counts * 10**6 / sampleDict[sample]
                geneObject.relativeExpression[sample] = relativeExpression

        # 统计每个exon的counts
        for geneObject in self.geneDict.values():
            for transcriptObject in geneObject.transcriptDict.values():
                for sample, counts in transcriptObject.countsExpression.items():
                    exonList = transcriptObject.exonList
                    for exon in exonList:
                        exonObject = self.exonDict[exon]
                        temp = exonObject.countsExpression.get(sample, 0)
                        temp = temp + counts
                        exonObject.countsExpression[sample] = temp
        # 计算每个exon的relativeExpression
        for exonObject in self.exonDict.values():
            for sample, counts in exonObject.countsExpression.items():
                relativeExpression = counts * 10**6 / sampleDict[sample]
                exonObject.relativeExpression[sample] = relativeExpression

        return None

    def computeCellLineExpression(self):
        '''
        change:
            计算每个gene, transcript, exon在每个细胞系中的相对表达量
                0. 处理gene
                1. 处理transcript
                    transcript在细胞系中的相对表达量 = 在一个细胞系的所有样本中的相对表达量的加和 / 细胞系中样本的数量
                    换句话说, transcript在细胞系中的相对表达量 = 细胞系中所有样本的平均相对表达量
                2. 处理exon
        '''

        # 处理gene
        for geneObject in self.geneDict.values():
            for cellLine, sampleList in self.celllineInfo.items():
                geneObject.cellLineExpression[cellLine] = numpy.mean([geneObject.relativeExpression.get(sample, 0) for sample in sampleList])

        # 处理transcript
        for geneObject in self.geneDict.values():
            for transcriptObject in geneObject.transcriptDict.values():
                for cellLine, sampleList in self.celllineInfo.items():
                    expression = numpy.mean([transcriptObject.relativeExpression.get(sample, 0) for sample in sampleList])
                    transcriptObject.cellLineExpression[cellLine] = expression

        # 处理exon
        for exonObject in self.exonDict.values():
            for cellLine, sampleList in self.celllineInfo.items():
                exonObject.cellLineExpression[cellLine] = numpy.mean([exonObject.relativeExpression.get(sample, 0) for sample in sampleList])

        return None

# 定义Exon
class Exon(object):
    '''
    属性
        基础属性
            status, str
            exonId, str
            exonVersion, int
            chr, str
            strand, str
            start, int, (与strand无关, start<end)
            end, int, (与strand无关, start<end)
            TSS, int, 转录起始位点
            TES, int, 转录终止位点
        额外属性
            countsExpression, dict, {<sample_name1>: <counts>, <sample_name2>: <counts>, ...}
            relativeExpression, dict, {<sample_name1>: <expression>, <sample_name2>: <expression>, ...}
            cellLineExpression, dict, {<cellLine1>: <expression>, <cellLine2>: <expression>, ...}
    方法
        _addTSSTES()                添加exon的TSS及TES信息
        dictGet()                   获取一个dict, key:value分别为exonId,chr,start,end,status
    '''
    def __init__(self, status, exonId, chr, strand, start, end, TSS=None, TES=None, exonVersion=-1):
        # 初始化属性
        self.status = status
        self.exonId = exonId
        self.exonVersion = int(exonVersion)
        self.chr = chr
        self.strand = strand
        self.start = start
        self.end = end
        self.TSS = TSS
        self.TES = TES
        # 额外属性
        self.countsExpression = {}
        self.relativeExpression = {}
        self.cellLineExpression = {}

    def _addTSSTES(self):
        '''
        change:
            判断该exon的TSS及TES
        '''
        if self.strand == '+':
            self.TSS = self.start
            self.TES = self.end
        else:
            self.TSS = self.end
            self.TES = self.start
        
        return None

    def dictGet(self):
        '''
        input:
            exonId, str
            chr, str
            strand, str
            start, int
            end, int
            status, str
        change:
            获取一个dict, key:value分别为exonId,chr,start,end,status
        output:
            dict
        '''
        temp = {"exonId": self.exonId,
                "chr": self.chr,
                "strand": self.strand,
                "start": self.start,
                "end": self.end,
                "status": self.status,
                "TSS": self.TSS,
                "TES": self.TES}

        return temp


# 定义Transcript
class Transcript(object):
    '''
    属性
        基础属性
            status, str
            transcriptId, str
            transcriptVersion, int
            transcriptName, str
            transcriptBiotype, str
            start, int, (与strand无关, start<end)
            end, int, (与strand无关, start<end)
            exonList, list, 该转录本的exon组成
            countsExpression, dict, {<sample_name1>: <counts>, <sample_name2>: <counts>, ...}
            relativeExpression, dict, {<sample_name1>: <expression>, <sample_name2>: <expression>, ...}
            cellLineExpression, dict, {<cellLine1>: <expression>, <cellLine2>: <expression>, ...}
        额外属性
            TSS, int, 转录起始位点
            TES, int, 转录终止位点
    方法
        _addTSSTES()        添加transcript的TSS及TES信息
        dictGet()           获取一个dict, 包含了该transcript的属性
        refresh()           重新整理transcript的exonList
    '''
    def __init__(self,status, transcriptId, transcriptName, transcriptBiotype, start, end, TSS=None, TES=None, transcriptVersion=-1, exonList=None, countsExpression=None, relativeExpression=None, cellLineExpression=None):
        self.status = status
        self.transcriptId = transcriptId
        self.transcriptVersion = int(transcriptVersion)
        self.transcriptName = transcriptName
        self.transcriptBiotype = transcriptBiotype
        self.start = start
        self.end = end
        self.exonList = (exonList, [])[exonList is None]
        self.countsExpression = (countsExpression, {})[countsExpression is None]
        self.relativeExpression = (relativeExpression, {})[relativeExpression is None]
        self.cellLineExpression = (cellLineExpression, {})[cellLineExpression is None]
        self.TSS = TSS
        self.TES = TES

    def _addTSSTES(self, strand):
        '''
        change:
            根据strand, 判断该transcript的TSS及TES
        '''
        strand = strand
        if strand == '+':
            self.TSS = self.start
            self.TES = self.end
        else:
            self.TSS = self.end
            self.TES = self.start

        return None

    def dictGet(self):
        '''
        change:
            获取一个dict, 包含了该transcript的属性
        output:
            dict, ...
        '''
        temp = {"transcriptId": self.transcriptId,
                "transcriptVersion": self.transcriptVersion,
                "transcriptName": self.transcriptName,
                "transcriptBiotype": self.transcriptBiotype,
                "status": self.status,
                "start": self.start,
                "end": self.end,
                "TSS": self.TSS,
                "TES": self.TES,
                "exonList": ", ".join(set(self.exonList))}

        return temp

    def refresh(self, exonExisted):
        '''
        input:
            exonExisted, dict, 参考Total对象中的exonExisted
        change:
            在total.exonDict及total.exonExisted被更改后
            可能该transcript的exonCombination会出现错误指向
            所以需要根据exonExisted更新exonCombination
                0. 遍历exonList, 寻找可进行映射的exon并进行映射
                1. 将exonList转为set格式后再转为list格式以去除重复
        return:
            bool, TRUE--exonCombination发生改变  False--exonCombination未发生改变
        '''
        exonExisted = exonExisted

        changeMarker = False
        # 对exon进行映射
        newExonList = self.exonList.copy()
        for i in range(0, len(newExonList)):
            exonId = newExonList[i]
            if exonId in exonExisted.keys():
                # 该exon可进行映射
                newExonList[i] = exonExisted[exonId]
                changeMarker = True

        # 去除exon重复
        newExonList = list(set(newExonList))
        self.exonList = newExonList

        if changeMarker is True:
            return True
        else:
            return False


# 定义Gene
class Gene(object):
    '''
    属性
        基础属性
            status, str,
            chr, str, gene所在的染色体号
            strand, str, gene的链的方向
            start, int, gene的start的位置(与strand无关, start<end)
            end, int, gene的end的位置(与strand无关, start<end)
            geneId, str
            geneVersion, int
            geneName, str
            geneBiotype, str
            exonList, list, 存储该gene所包含的所有exonId
            transcriptDict, dict, 存储映射的Transcript对象, key为id
            transcriptIndex, dict, 建立基因中transcript的索引 {<start>: {<end>: <transcriptId>}}
            transcriptExisted, dict, 存储transcriptId的映射关系
        额外属性
            TSS, int, 转录起始位点
            TSE, int, 转录终止位点
            geneClass, str, TSS-PAS or ATSS-PAS or TSS-APA or ATSS-APA
            countsExpression, dict, {<sample_name1>: <counts>, <sample_name2>: <counts>, ...}
            relativeExpression, dict, {<sample_name1>: <expression>, <sample_name2>: <expression>, ...}
            cellLineExpression, dict, {<cellLine1>: <expression>, <cellLine2>: <expression>, ...}
    方法
        _exonCheck()            检查exon是否已被记录
        _transcriptCheck()      检查transcript是否已被记录
        _transcriptExistedAdd() 向transcriptExisted中添加映射
        _transcriptMerge()      合并同一gene中的两个transcript
        _addTSSTES()            添加gene的TSS及TES信息
        dictGet()               获取一个dict, 含有该gene的基本属性
        exonUseGet()            查询该gene中使用了指定的exon的transcript
        transcriptAdd()         向transcriptDict中添加Transcript对象
        transcriptIdGet()       输入一个transcriptId或transcriptName, 返回一个考虑映射关系后的transcriptId
        transcriptIndexBuild()  建立符合过滤条件的transcript的index
        refresh()               更新gene的exonList, 并更新transcriptDict中每个transcript的exonList
        reIndex()               重新建立transcriptIndex
    '''
    def __init__(self, status, chr, strand, start, end, geneId, geneName, geneBiotype, TSS=None, TES=None, geneVersion=-1, exonList=None, transcriptDict=None, transcriptIndex=None, transcriptExisted=None):
        # 基础属性
        self.status = status
        self.chr = chr
        self.strand = strand
        self.start = start
        self.end = end
        self.geneId = geneId
        self.geneVersion = int(geneVersion)
        self.geneName = geneName
        self.geneBiotype = geneBiotype
        self.exonList = (exonList, [])[exonList is None]
        self.transcriptDict = (transcriptDict, {})[transcriptDict is None]
        self.transcriptIndex = (transcriptIndex, {})[transcriptIndex is None]  # {<start>: {<end>: [], ...}, ...}
        self.transcriptExisted = (transcriptExisted, {})[transcriptExisted is None]  # {<old transcript id>: <existed transcript id>}
        self.TSS = TSS
        self.TES = TES
        # 额外属性
        self.geneClass = None  # TSS-PAS or ATSS-PAS or TSS-APA or ATSS-APA
        self.countsExpression = {}
        self.relativeExpression = {}
        self.cellLineExpression = {}

    def _exonCheck(self, exonId):
        '''
        input:
            exonId, str
        change:
            检查该gene中是否已记录指定的exonId
        output:
            bool, True-已记录 False-未记录
        '''
        exonId = exonId
        
        if exonId in self.exonList:
            return True
        else:
            return False

    def _transcriptCheck(self, status, transcriptId, start, end, exonCombination):
        '''
        input:
        change:
            检查transcript是否已被记录
        output:
            list, transcriptId已映射--[0, <existed transcript id>]--旧
                  transcriptId已记录--[1, None]--旧
                  transcriptId未记录 且start未记录--[2, None]--新
                  transcriptId未记录 且start已记录 且end未记录--[3, None]--新
                  transcriptId未记录 且start已记录 且end已记录 且exon组成与其他transcript一致--[4, <existed transcript id>]--旧
                  transcriptId未记录 且start已记录 且end已记录 且exon组成独一无二--[5, None]--新

        '''
        status = status
        transcriptId = transcriptId
        start = start
        end = end
        exonCombination = set(exonCombination)

        # 检查该transcriptId是否映射到了其他transcriptId
        if transcriptId in self.transcriptExisted.keys():
            # 已映射, 认为该transcript已存在
            existedTranscriptId = self.transcriptExisted[transcriptId]
            return [0, existedTranscriptId]
        else:
            # 未映射到其他transcriptId
            pass
        # 检查在transcript对应的gene中该transcriptId是否已存在
        if transcriptId in self.transcriptDict.keys():
            # 该transcriptId已记录
            return [1, None]
        else:
            # 该transcriptId未记录
            # 检查transcript的start是否已记录
            if start not in self.transcriptIndex.keys():
                # transcriptId未记录 且start未记录
                return [2, None]
            else:
                # transcriptId未记录 且start已记录
                # 检查transcript的end是否已记录
                if end not in self.transcriptIndex[start]:
                    # transcriptId未记录 且start已记录 且end未记录
                    return [3, None]
                else:
                    # transcriptId已记录 且start已记录 且end已记录
                    # 检查transcript的exon组成是否与相应的transcript相同
                    existedTranscriptIdList = self.transcriptIndex[start][end]
                    for existedTranscriptId in existedTranscriptIdList:
                        existedExonCombination = self.transcriptDict[existedTranscriptId].exonList
                        existedExonCombination = set(existedExonCombination)
                        # 存在transcript的exon组成与该transcript一致
                        if exonCombination == existedExonCombination:
                            return [4, existedTranscriptId]
                    # 没有transcript的exon组成与该transcript一致
                    return [5, None]

    def _transcriptExistedAdd(self, oldTranscriptId, newTranscriptId):
        '''
        change:
            向transcriptExisted中添加映射{oldTranscriptId: newTranscriptId}
        return:
            bool, True-已添加映射  False-自我映射
        '''
        oldTranscriptId = oldTranscriptId
        newTranscriptId = newTranscriptId

        keys = set(self.transcriptExisted.keys())
        values = set(self.transcriptExisted.values())

        # 寻找到所有可以映射到newTranscriptId的transcriptId
        oldTranscriptIdSet = {oldTranscriptId}
        workSet = {oldTranscriptId}
        while workSet:
            transcriptId = workSet.pop()
            if transcriptId in values:
                # 有其他的transcriptId映射到了该transcriptId
                temp = {k for k, v in self.transcriptExisted.items() if v==transcriptId}
                # workSet添加这些transcriptId, 检查是否还有transcriptId可映射到这些transcriptId        
                workSet = workSet.union(temp)
            else:
                oldTranscriptIdSet.add(transcriptId)
        # 寻找到最终的newTranscriptId
        while newTranscriptId in keys:
            newTranscriptId = self.transcriptExisted[newTranscriptId]

        # 不进行自我映射
        if oldTranscriptId == newTranscriptId:
            return False
        elif newTranscriptId in oldTranscriptIdSet:
            return False

        for transcriptId in oldTranscriptIdSet:
            self.transcriptExisted[transcriptId] = newTranscriptId

        return True

    def _transcriptMerge(self, transcriptId1, transcriptId2):
        '''
        change:
            合并两transcript
            0. 检查两transcript是否存在, 若存在则进行下一步
            1. 确定映射方向并映射
                若2个transcript的status为KNOWN, 则在self.transcriptExisted中将transcriptId2映射到transcriptId1
                若1个transcript的status为KNOWN, 则在self.transcriptExisted中将非KNOWN的transcriptId映射到KNOWN的transcriptId
                若0个transcript的status为KNOWN, 则在self.transcriptExisted中将transcriptId2映射到transcriptId1
            2. 将deletedTranscript的countsExpression添加到remainedTranscript的countsExpression中
            3. 在self.transcriptDict中删除deletedTranscript
            4. 重建transcriptIndex
            # 5. 不更新gene的exonList
        return:
            tuple, 正常情况下会返回(<deletedTranscriptId>, <remainedTranscriptId>), 但当deletedTranscriptId==remainedTranscriptId时会返回("", "")
        '''
        transcriptId1 = transcriptId1
        transcriptId2 = transcriptId2

        deletedTranscriptId = None
        remainedTranscriptId = None

        # 检查两transcript是否存在
        if transcriptId1 not in self.transcriptDict.keys():
            raise ValueError("[Error]_transcriptMerge() the transcriptId1 {} is not recorded".format(transcriptId1))
        elif transcriptId2 not in self.transcriptDict.keys():
            raise ValueError("[Error]_transcriptMerge() the transcriptId2 {} is not recorded".format(transcriptId2))
        else:
            pass

        # 确定映射方向
        transcriptObject1 = self.transcriptDict[transcriptId1]
        transcriptObject2 = self.transcriptDict[transcriptId2]
        if transcriptObject1.status=="KNOWN" and transcriptObject2.status=="KNOWN":
            # 将transcriptId2映射到transcriptId1
            deletedTranscriptId = transcriptId2
            remainedTranscriptId = transcriptId1
        elif transcriptObject1.status=="KNOWN":
            # 将transcriptId2映射到transcriptId1
            deletedTranscriptId = transcriptId2
            remainedTranscriptId = transcriptId1
        elif transcriptObject2.status=="KNOWN":
            # 将transcriptId1映射到transcriptId2
            deletedTranscriptId = transcriptId1
            remainedTranscriptId = transcriptId2
        else:
            # 将transcriptId2映射到transcriptId1
            deletedTranscriptId = transcriptId2
            remainedTranscriptId = transcriptId1
        marker = self._transcriptExistedAdd(oldTranscriptId=deletedTranscriptId, newTranscriptId=remainedTranscriptId)

        # 处理异常: deletedTranscriptId == remainedTranscriptId
        if marker is None:
            return ("", "")

        deletedTranscript = self.transcriptDict[deletedTranscriptId]

        # 更新countsExpression
        for sample, counts in deletedTranscript.countsExpression.items():
            temp = self.transcriptDict[remainedTranscriptId].countsExpression.get(sample, 0)
            temp = temp + counts
            self.transcriptDict[remainedTranscriptId].countsExpression[sample] = temp

        # 修改transcriptIndex
        try:
            self.transcriptIndex[self.transcriptDict[deletedTranscriptId].start][self.transcriptDict[deletedTranscriptId].end] = remainedTranscriptId
        except KeyError as k:
            # 目前这个地方出现问题的原因:
            # 在reIndex()时, transcript的start/end与index中的start/end不一致(校正后改变了start/end的位置)
            k
            pass

        # 删除deletedTranscript
        self.transcriptDict.pop(deletedTranscriptId)

        return (deletedTranscriptId, remainedTranscriptId)

    def _addTSSTES(self):
        '''
        change:
            根据strand, 判断该gene及所包含的transcript的TSS及TES
        '''
        if self.strand == '+':
            self.TSS = self.start
            self.TES = self.end
        else:
            self.TSS = self.end
            self.TES = self.start

        for transcriptId, transcriptObject in self.transcriptDict.items():
            transcriptObject._addTSSTES(strand=self.strand)
        
        return None

    def dictGet(self):
        '''
        change:
            获取一个dict, 含有该gene的基本属性
        output:
            dict, ...
        '''

        temp = {"geneId": self.geneId,
                "geneVersion": self.geneVersion,
                "geneName": self.geneName,
                "geneBiotype": self.geneBiotype,
                "status": self.status,
                "strand": self.strand,
                "chr": self.chr,
                "start": self.start,
                "end": self.end,
                "TSS": self.TSS,
                "TES": self.TES,
                "transcript": ", ".join(set(self.transcriptDict.keys())),
                "exon": ", ".join(set(self.exonList))}

        return temp

    def exonUseGet(self, exonId):
        '''
        change:
            查询该gene中使用了指定的exon的transcript
        return:
            dict, 若该gene中仅含有指定exon而无transcript使用该exon, 则返回{<geneId>: set()}
                  若该gene中存在transcript使用了该exon, 则返回{<geneId>: {<transcriptId>, ...}}
                  若该gene中不存在指定exon, 则返回{}
        '''
        exonId = exonId

        transcriptSet = set()
        for transcriptId, transcriptObject in self.transcriptDict.items():
            if exonId in transcriptObject.exonList:
                transcriptSet.add(transcriptId)

        if len(transcriptSet) == 0:
            # 没有transcript使用该exon
            if exonId in self.exonList:
                # 该gene中仅含有指定exon而无transcript使用该exon
                return {self.geneId: set()}
            else:
                # 该gene中不存在指定exon
                return {}
        else:
            # 存在transcript使用该exon
            return {self.geneId: transcriptSet}

    def transcriptAdd(self, status, transcriptId, transcriptVersion, transcriptName, transcriptBiotype, start, end, exonList, sample, counts):
        '''
        input:
            status, str
            transcriptId, str
            transcriptName, str
            transcriptBiotype, str
            start, int, (与strand无关, start<end)
            end, int, (与strand无关, start<end)
            exonList, list, 该转录本的exon组成
            sample, str, 当前transcript所在的样本的样本名
            counts, int, 在quantification文件中该transcript的counts数
        change:
            添加transcript
        output:
            bool, None
        '''
        status = status
        transcriptId = transcriptId
        transcriptName = transcriptName
        transcriptBiotype = transcriptBiotype
        start = start
        end = end
        exonList = exonList
        sample = sample
        counts = int(counts)

        checkMarker = self._transcriptCheck(status=status, transcriptId=transcriptId, start=start, end=end, exonCombination=exonList)
        if checkMarker[0] == 0:
            # 已映射 旧
            existedTranscriptId = checkMarker[1]
            temp = self.transcriptDict[existedTranscriptId].countsExpression.get(sample, 0)
            self.transcriptDict[existedTranscriptId].countsExpression[sample] = counts + temp
            return None
        elif checkMarker[0] == 1:
            # transcriptId已记录 旧
            self.transcriptDict[transcriptId].countsExpression[sample] = counts
            return None
        elif checkMarker[0] == 2:
            # transcriptId未记录 且start未记录 新
            self.transcriptDict[transcriptId] = Transcript(status=status,
                                                           transcriptId=transcriptId,
                                                           transcriptVersion=transcriptVersion,
                                                           transcriptName=transcriptName,
                                                           transcriptBiotype=transcriptBiotype,
                                                           start=start,
                                                           end=end,
                                                           exonList=exonList)
            self.transcriptDict[transcriptId].countsExpression[sample] = counts
            self.transcriptIndex[start] = {end: [transcriptId]}
            return None
        elif checkMarker[0] == 3:
            # transcriptId未记录 且start已记录 且end未记录 新
            self.transcriptDict[transcriptId] = Transcript(status=status,
                                                           transcriptId=transcriptId,
                                                           transcriptVersion=transcriptVersion,
                                                           transcriptName=transcriptName,
                                                           transcriptBiotype=transcriptBiotype,
                                                           start=start,
                                                           end=end,
                                                           exonList=exonList)
            self.transcriptDict[transcriptId].countsExpression[sample] = counts
            self.transcriptIndex[start][end] = [transcriptId]
            return None
        elif checkMarker[0] == 4:
            # transcriptId未记录 且start已记录 且end已记录 且exon组成与其他transcript一致 旧
            existedTranscriptId = checkMarker[1]
            temp = self.transcriptDict[existedTranscriptId].countsExpression.get(sample, 0)
            self.transcriptDict[existedTranscriptId].countsExpression[sample] = counts + temp
            self._transcriptExistedAdd(oldTranscriptId=transcriptId, newTranscriptId=existedTranscriptId)
            return None
        elif checkMarker[0] == 5:
            # transcriptId未记录 且start已记录 且end已记录 且exon组成独一无二 新
            self.transcriptDict[transcriptId] = Transcript(status=status,
                                                           transcriptId=transcriptId,
                                                           transcriptVersion=transcriptVersion,
                                                           transcriptName=transcriptName,
                                                           transcriptBiotype=transcriptBiotype,
                                                           start=start,
                                                           end=end,
                                                           exonList=exonList)
            self.transcriptDict[transcriptId].countsExpression[sample] = counts
            self.transcriptIndex[start][end].append(transcriptId)
            return None
        else:
            raise ValueError("[Error]transcriptAdd() the checkMarker is {}, should be 0,1,2,3,4,5".format(checkMarker))

    def transcriptIdGet(self, transcriptId=None, transcriptName=None):
        '''
        input:
            transcriptId, str, =None
            transcriptName, str, =None
        change:
            获取gene中的transcript id(已考虑transcriptExisted中的映射)
        output:
            str or bool, 若能寻找到transcriptId则返回transcriptId, 否则返回None
        '''
        transcriptId = transcriptId
        transcriptName = transcriptName

        if transcriptId is None and transcriptName is None:
            raise ValueError("[Error]transcriptIdGet() need transcriptId or transcriptName")
        elif transcriptName is None:
            # 根据transcriptId进行寻找
            if transcriptId in self.transcriptExisted.keys():
                # transcriptId被映射到了其他transcriptId
                transcriptId = self.transcriptExisted[transcriptId]
                return transcriptId
            else:
                # transcriptId未映射到其他transcriptId
                if transcriptId in self.transcriptDict.keys():
                    # 该transcriptId已被记录
                    return transcriptId
                else:
                    # 该transcriptId未被记录
                    return None
        else:
            # 根据transcriptName进行寻找
            for transcriptId, transcriptObject in self.transcriptDict.items():
                if transcriptName == transcriptObject.transcriptName:
                    # 已根据transcriptName寻找到对应的transcriptId
                    return transcriptId
            # 未能根据transcriptName寻找到transcriptId
            return None

    def transcriptIndexBuild(self, siteType="start", statusSet=None):
        '''
        input:
            siteType, str, "start" or "end" or "TSS" or "TES"
            statusSet, set, 只有status在这个set中的transcript才会被纳入index
        return:
            index, dict, {<siteType>: {siteType2: geneId}, ...}
        '''
        siteType = siteType
        statusSet = (statusSet, ("KNOWN", "NOVEL", "CORRECTED-start", "CORRECTED-end", "CORRECTED-start-end", "CORRECTED-end-start"))[statusSet is None]

        mapped = {"start": "end", "end": "start", "TSS": "TES", "TES": "TSS"}
        index = {}
        for transcriptId, transcriptObject in self.transcriptDict.items():
            status = transcriptObject.status
            # 过滤掉不符合条件的gene
            if status not in statusSet:
                continue
            data = {"start": transcriptObject.start, "end": transcriptObject.end, "TSS": transcriptObject.TSS, "TES": transcriptObject.TES}
            if data[siteType] not in index.keys():
                index[data[siteType]] = {data[mapped[siteType]]: [transcriptId]}
            else:
                if data[mapped[siteType]] not in index[data[siteType]].keys():
                    index[data[siteType]][data[mapped[siteType]]] = [transcriptId]
                else:
                    index[data[siteType]][data[mapped[siteType]]].append(transcriptId)

        return index

    def refresh(self, exonExisted):
        '''
        input:
            exonExisted, dict, 参考Total对象中的属性exonExisted
        change:
            更新gene的exonList, 并更新transcriptDict中每个transcript的exonList
                0. 更新gene中每个transcript的exonList
                1. 更新gene的exonList
                    将每个transcript的exon都添加到gene的exonList中
                    去除gene的exonList中的重复项
                    将gene的exonList中可映射的exon进行映射
        return:
            list, True--该gene的exonList已更改 False--该gene的exonList未更改  [<bool>, <changed transcriptId>]
        '''
        exonExisted = exonExisted

        # 更新gene中每个transcript的exonList
        changedTranscriptId = []  # 记录exonList已更改的transcriptId
        for transcriptId, transcriptObject in self.transcriptDict.items():
            temp = transcriptObject.refresh(exonExisted=exonExisted)
            if temp is True:
                # 该transcript的exonList已更改
                changedTranscriptId.append(transcriptId)

        # 将每个transcript的exon都添加到gene的exonList中
        newExonList = self.exonList.copy()
        for transcriptId, transcriptObject in self.transcriptDict.items():
            newExonList = newExonList + transcriptObject.exonList
        # 去除gene的exonList中的重复项
        newExonList = list(set(newExonList))

        # 将gene的exonList中可映射的exon进行映射
        changeMarker = False
        for i in range(0, len(newExonList)):
            exonId = newExonList[i]
            if exonId in exonExisted.keys():
                # 该exon可进行映射
                newExonList[i] = exonExisted[exonId]
        if set(newExonList) != set(self.exonList):
            changeMarker = True

        # 更新gene的exonList
        if changeMarker is True:
            self.exonList = list(set(newExonList))
            return [True, changedTranscriptId]
        else:
            return [False, changedTranscriptId]

    def reIndex(self):
        '''
        change:
            遍历self.transcriptDict, 重新建立self.transcriptIndex
                0. 声明dict类型变量newIndex以存储新的index
                1. 遍历每一个transcript对象
                2. 检查transcript的start是否已记录
                    若start未记录, 则在newIndex添加该transcript, 结束函数
                    否则进行下一步
                3. 检查transcript的end是否已记录
                    若end未记录, 则在newIndex添加该transcript, 结束函数
                    否则进行下一步
                4. 检查transcript的exonCombination是否已记录
                    若exonCombination未记录, 则在newIndex添加该transcript, 结束函数
                    否则表明transcript与newIndex中的transcript存在重复, 合并transcript, 结束函数
        return:
            None
        '''

        # 存储新的index
        newIndex = {}
        # 遍历每一个transcript对象
        tempTranscriptDict = self.transcriptDict.copy()
        for transcriptId, transcriptObject in tempTranscriptDict.items():
            breakMarker = False  # 表示该transcript已判断完毕
            start = transcriptObject.start
            end = transcriptObject.end
            exonCombination = set(transcriptObject.exonList)
            # 检查transcript的start是否已记录
            if start not in newIndex.keys():
                # start未记录
                newIndex[start] = {end: [transcriptId]}
            else:
                # start已记录
                # 检查transcript的end是否已记录
                if end not in newIndex[start].keys():
                    # end未记录
                    newIndex[start][end] = [transcriptId]
                else:
                    # end已记录
                    # 检查transcript的exonCombination是否已记录
                    tempTranscriptList = newIndex[start][end]
                    for id in tempTranscriptList:
                        if id in self.transcriptExisted.keys():
                            id = self.transcriptExisted[id]
                        tempExonCombination = set(self.transcriptDict[id].exonList)
                        if exonCombination == tempExonCombination:
                            # transcript存在重复
                            self._transcriptMerge(transcriptId1=id, transcriptId2=transcriptId)
                            breakMarker = True
                            break
                    if breakMarker is True:
                        continue
                    else:
                        # exonCombination未记录
                        temp = newIndex[start][end]
                        temp.append(transcriptId)
                        newIndex[start][end] = temp
        # 保存新的Index
        self.transcriptIndex = newIndex

        return None
