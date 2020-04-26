from itertools import combinations
from Node import Cell
from anytree import Node, RenderTree
import re
i = 0

global p
global parent
global alignmentTree
global tempAlignmentTree
global findingPaths# backtrack ile bulunan yollar
global gapScore
global isMultipleAlignment
global indexofProfile # meanin icini doldurmak icin index olarak kullandım
global resultOfMultipleAlignment
resultOfMultipleAlignment= dict()
global indexForResultofMultipleAlignment ## multiple alingment resultunu doldurmak için
isMultipleAlignment=0

gapScore = input("Please, enter gap score:")

if(gapScore.lstrip("-").isdigit()):
    gapScore = int(gapScore)
else:
    gapScore = input("Please, enter invalid gap score:")
    while(gapScore.lstrip("-").isdigit() == 0):
        gapScore = input("Please, enter invalid gap score:")
    gapScore = int(gapScore)

parent = 0
alignmentTree = dict()
tempAlignmentTree = dict()
scoreList = list()
directionList = list()
mean = dict()
firstSequence = ''
SecondSequence = ''

#read input file
file_object  = open("..\Input1.txt", "r")
fileArr = file_object.read()

#read blosum62
blosum62  = open("..\Blosum62.txt", "r")
blosumArr = blosum62.read()
temp = ""
Blosum62List = list()
Blosum62List = blosumArr.split('\n')

file_sequence = list()
file_sequence = fileArr.split("\n")
sequences = list()
sequenceName = list()

for j in file_sequence:
    if (('>' not in j) and (not j.__eq__(''))):
        sequences.append(j)
    elif('>' in j):
        sequenceName.append(j)

print("SEQUENCES")
print("-----------------")
for a in range(0, len(sequences)):
    print(sequenceName[a],sequences[a])

combOfSeq = combinations(sequences, 2)

for x in range(0, Blosum62List.__len__()):
    for char in Blosum62List[x]:
        if char != " ":
            temp += char
    Blosum62List[x] = temp
    temp = ""

def getScoreFromBlosum62(firstSequence, SecondSequence):
    m = 0
    countFirst = 0
    countSecond = 0

    for first in Blosum62List[0]:
        if(first == firstSequence):
            countFirst = m
            break
        else:
            m += 1
    m = 1
    for x in range(1, Blosum62List.__len__()):
        if(SecondSequence == Blosum62List[x][0]):
            countSecond = m
            break
        else:
            m += 1
    m = 1
    x = 1
    while x <= countFirst:
        if(Blosum62List[countSecond][m] == "-"):
            m += 1
        else:
            x += 1
            m += 1

    if(Blosum62List[countSecond][m] == "-"):
        return int(Blosum62List[countSecond][m+1])*(-1)
    else:
        return int(Blosum62List[countSecond][m])

def reverse(s):
  str = ""
  for i in s:
    str = i + str
  return str

#Exact match oranı buluyo
global maxScore
maxScore=list()
global alignSequence
alignSequence = list()
def alignMostSimilarSequences(sequence1,sequence2):
    alignSequence.clear()
    """en yuksek oranda benzerlik iceren sequence'lar icin birlestirme yapiyor"""
    # gelen sequence'larin list tipine donusturulmesi gerekiyor !!
    seq1=list(sequence1)
    seq2=list(sequence2)
    tempSeq = ""
    for i in range(len(seq1)):
        if (len(seq1[i]) == 1 and len(seq2[i]) == 1):  # sadece 1 tane aa varsa
            if (seq2[i].__eq__(seq1[i])):  # aminoacidler birb esit ise direk bir tanesini ekliyoruz
                alignSequence.append("*"+seq1[i]+"*")
            else:  # aminoacidler birb esit degilse ikisinin arasina virgul koyarak ekliyoruz
                temp = seq1[i] + "," + seq2[i]
                alignSequence.append("*"+temp+"*")
        else:  # eger eslesecek olan aminoacd biden fzla ise
            if (seq2[i].__eq__(seq1[i])):  # ilk basta ayni mi diye bakiyoruz
                alignSequence.append("*"+seq1[i]+"*")  # ayni ise direk ekliyoruz
            else:  # ayni degilse for ile donmek zorundayiz
                if (len(seq1[i]) > len(seq2[i])):
                    for sLen1 in range(len(seq1[i])):
                        for sLen2 in range(len(seq2[i])):
                            if (seq1[i][sLen1].__eq__(seq2[i][sLen2])):
                                tempSeq = tempSeq + seq1[i][sLen1]
                                break
                            elif (sLen2 == len(seq2[i]) - 1):
                                tempSeq = tempSeq + seq1[i][sLen1]
                    for sLen2 in range(len(seq2[i])):
                        if (seq2[i][sLen2] not in tempSeq):
                            tempSeq = tempSeq + "," + seq2[i][sLen2]
                    alignSequence.append("*"+tempSeq+"*")
                elif (len(seq1[i]) < len(seq2[i]) or len(seq1[i]) == len(seq2[i])):
                    for sLen2 in range(len(seq2[i])):
                        for sLen1 in range(len(seq1[i])):
                            if (seq1[i][sLen1].__eq__(seq2[i][sLen2])):
                                tempSeq = tempSeq + seq1[i][sLen1]
                                break
                            elif (sLen1 == len(seq1[i]) - 1):
                                tempSeq = tempSeq + seq2[i][sLen2]
                    for sLen1 in range(len(seq1[i])):
                        if (seq1[i][sLen1] not in tempSeq):
                            tempSeq = tempSeq + "," + seq1[i][sLen1]
                    alignSequence.append("*"+tempSeq+"*")
    return alignSequence

def multipleAlignment():
    global isMultipleAlignment
    global indexForResultofMultipleAlignment
    global indexofProfile
    global resultOfMultipleAlignment
    x=""
    #tum seq sequences icinde tutuluyor(liste seklinde)
    #sequenceValue = dict() seq larin score ve align sekilleri
    isMultipleAlignment=1# multiple alignmnet kosuluna geldiysek true
    global alignmentTree
    for i in alignmentTree:
        indexofProfile = i
        score = 0
        s1 = ""
        s2 = ""
        #mean dict olarak tuttugumuz degisken hangi sayinin hangi seq a denk geldigini tutuyor
        value = alignmentTree[i]
        #birbirine align edilecek sequence lari dictionaryden alip first ve second sequence a yazdiriyoruz
        value = value[1:-1]# bastaki ve sinrdaki parantezleri silmek icin
        firstSequence,secondSequence=value.split(",")
        if(firstSequence[0].__eq__("s")):# seq ile seq align ediliyorsa
            a=int(firstSequence[1])-1
            firstSequence=sequences[a]
            b = int(secondSequence[1]) - 1
            secondSequence = sequences[b]
            for key in sequenceValue:
                temp = sequenceValue[key].split(",")
                if((temp[0].__eq__(firstSequence) and temp[1].__eq__(secondSequence)) or (temp[1].__eq__(firstSequence) and temp[0].__eq__(secondSequence))):
                    if(score<key):
                        score=key
                        s1=temp[2]
                        s2=temp[3]
            new=alignMostSimilarSequences(s1, s2)
            bag="s"+str(a+1)
            resultOfMultipleAlignment[bag]= s1
            bag = "s" + str(b + 1)
            resultOfMultipleAlignment[bag] = s2
            x=''.join(new)
            v=x.split("*")
            last=list()
            for q in v:
                if(q is not ''):
                    last.append(q)
            mean[i] =last
        elif(not firstSequence[0].__eq__("s")):# seq ile profil yada profil ile profil align ediliyorsa
            profileFirst = mean[int(firstSequence[0])]# ilk profilim
            if(secondSequence[0].__eq__("s")): #profil ile seq align ediliyosa
                sequenceSecond = sequences[int(secondSequence[1])-1] # seq
                indexForResultofMultipleAlignment = secondSequence[1] # multiple alignment resultunu doldurmak için
                getAlignment(profileFirst,sequenceSecond)  #elde ettigim profile ve sekansı alingmenta gonderdim
            elif(not secondSequence[0].__eq__("s")): #profil profil align ediyosam
                profileSecond = mean[int(secondSequence[0])]  # ikinci profilim
                getAlignment(profileFirst, profileSecond)

global sequenceValue # hangi seq'larin hangi score ile alignment olduklarini tutuyor
sequenceValue = dict()

def getSimilarityScore(originalS1,originalS2,sequence1,sequence2):
        a=0
        b=sequence1.__len__()
        for x, y in zip(sequence1, sequence2):
             if x == y:
                 a+=1
        similarity = a/b
        maxScore.append(similarity)

        key=similarity
        value=originalS1+","+originalS2+","+sequence1+","+sequence2
        sequenceValue.update({key:value})

findingPaths=list()

def backtrack(nodes):
    # diagonal=0 top=1 left=2
    tempSequence=""# bulunan 1 yol
    tempNode=nodes# asil nodu kaybetmemek icin
    i = len(nodes)-1
    j = len(nodes[0])-1
    a=0# 2. dallanmanin oldugu yerde 1. dallanmanin i sini tutmak icin
    b=0
    dimention=0 # 2. dallanmanin oldugu yerde 1. dallanmanin yonunu tutmak icin
    while i >= 0 and j >= 0:
            num = 0  # kac tane dallanma oldugunu tutacak
            if(i == 0 and j == 0):
                break
            if(tempNode[i][j].get_top() == 1):
               num=num+1
            if(tempNode[i][j].get_left() == 1):
                num = num + 1
            if(tempNode[i][j].get_diagonal() == 1):
                num = num + 1
            if(num==1):
                if (tempNode[i][j].get_top() == 1):
                    tempSequence=tempSequence+"1"
                    i-=1
                elif (tempNode[i][j].get_left() == 1):
                    tempSequence = tempSequence + "2"
                    j-=1
                elif (tempNode[i][j].get_diagonal() == 1):
                    tempSequence = tempSequence + "0"
                    i-=1
                    j-=1

            elif (num==2 or num==3):
                if(a !=0 and b!=0):
                    if(dimention==0):
                        tempNode[a][b].set_diagonal(1)
                    elif (dimention == 1):
                        tempNode[a][b].set_top(1)
                    elif (dimention == 2):
                         tempNode[a][b].set_left(1)
                a = i
                b = j

                if (tempNode[i][j].get_top() == 1):
                    tempSequence=tempSequence+"1"
                    tempNode[i][j].set_top(0)
                    dimention=1
                    i-=1
                elif (tempNode[i][j].get_left() == 1):
                    tempSequence = tempSequence + "2"
                    tempNode[i][j].set_left(0)
                    j-=1
                    dimention = 2
                elif (tempNode[i][j].get_diagonal() == 1):
                    tempSequence = tempSequence + "0"
                    tempNode[i][j].set_diagonal(0)
                    dimention = 0
                    i-=1
                    j-=1

    if(isMultipleAlignment==1):
        findingPaths.clear()
    findingPaths.append(tempSequence)
    if (isMultipleAlignment==0):
        temp=findingPaths[-2:]
        if(len(temp)>1):
            if(temp[0].__eq__(temp[1])):
                del findingPaths[-1]
                return
        backtrack(tempNode)


def getMaxAlignmentScore(nodes,sequence1,sequence2):
    # backtrack ten bulunan birden fazla path icin maximum benzerlik oranini bulmaya yarayan fonk
    for path in findingPaths:
        s1 = ""
        s2 = ""
        i = len(nodes) - 1
        j = len(nodes[0]) - 1
        for k in range(len(path)):
            if(path[k].__eq__("0")): #diagonal
                s1=s1+sequence1[j-1]
                s2=s2+sequence2[i-1]
                i-=1
                j-=1
            elif(path[k].__eq__("1")): #top
                s1=s1+"-"
                s2=s2+sequence2[i-1]
                i-=1
            elif(path[k].__eq__("2")):# left
                s1=s1+sequence1[j-1]
                s2=s2+"-"
                j-=1
        s1=reverse(s1)
        s2=reverse(s2)
        getSimilarityScore(sequence1, sequence2, s1, s2)

def getProfileForMultipleAlignment(nodes,sequence1,sequence2): #backtrakten elde edilen yol ile profile bulmak için
    global indexForResultofMultipleAlignment
    global resultOfMultipleAlignment
    output = str(type(sequence2)).split("'") ## eger sequence2 profil değilse seq ise (str ise) onun align edilmiş halini saklamam lazım
    sequence2Type = output[1] # str ise sequence demektir (profile değil)

    s1 = list()
    s2 = list()
    i = len(nodes) - 1
    j = len(nodes[0]) - 1
    path = findingPaths[0] #multiple alignmentta backtrackten gelen ilk bulduğum yol
    for k in range(path.__len__()):
        if(path[k].__eq__("0")): #diagonal
            s1.append(sequence1[j-1])
            s2.append(sequence2[i-1])
            i-=1
            j-=1
        elif(path[k].__eq__("1")): #top
            s1.append("-")
            s2.append(sequence2[i-1])
            i-=1
        elif(path[k].__eq__("2")):# left
            s1.append(sequence1[j-1])
            s2.append("-")
            j-=1
    s1=list(reversed(s1))
    s2=list(reversed(s2))
    if (sequence2Type.__eq__("str")):
        seq2=''.join(s2)
        resultOfMultipleAlignment["s"+indexForResultofMultipleAlignment] = seq2

    profile = list()
    lenght=len(s1)
    for x in range(lenght):
        item1=s1[x].split(",")
        item2=s2[x].split(",")
        for i in range(len(item1)):
            for j in range(len(item2)):
                if(item1[i] == item2[j]):
                    item2[j]=""
        combination=""
        for k in range(len(item1)):
            if(item1[k]!=""):
                combination+=item1[k]
                combination+=","
        for l in range(len(item2)):
            if (item2[l] != ""):
                combination += item2[l]
                combination+=","
        combination=combination[:-1]
        profile.append(combination)
    mean[indexofProfile]=profile

    if(len(resultOfMultipleAlignment) == len(sequences)):
        print("MULTIPLE ALIGNMENT RESULT")
        print("-------------------------")
        for align in resultOfMultipleAlignment:
            temp = align.split("s")
            for name in range(0, len(sequenceName)):
                if((int(temp[1]) -1) == name):
                    print(sequenceName[name], resultOfMultipleAlignment[align])

def getAlignment(sequence1,sequence2):
    x = sequence1.__len__()+1
    y = sequence2.__len__()+1

    #ilki colon sayisi ikincisi satir sayisi
    nodes = [[Cell() for l in range(x)] for r in range(y)]
    for i in range(y):#sequence 2 icin satir satir bakmak
        for j in range(x):  # sequence 1 icin sutun sutun bakmak
            if i==0:#ilk satira bakiyoruz
                if j == 0:# ilk satir ve ilk sutunda bulunan eleman icin->0 olan
                    nodes[i][j].set_score(0)
                else:# ilk eleman degilse
                    if(isMultipleAlignment==1):
                        element = sequence1[j-1].split(",")
                        tempGapScore =element.__len__() * gapScore
                        leftScore = nodes[i][j - 1].get_score() + tempGapScore
                        nodes[i][j].set_score(leftScore)
                        nodes[i][j].set_left(1)
                    else:
                        leftScore=nodes[i][j].get_score()+gapScore
                        nodes[i][j].set_score(leftScore)
                        nodes[i][j].set_left(1)
            else:#ilk satir degilse
                if j == 0:#ilk sutuna bakiyorsak
                    if (isMultipleAlignment == 1):
                        element = sequence2[i-1].split(",")
                        tempGapScore = element.__len__() * gapScore
                        topScore = nodes[i-1][j].get_score() + tempGapScore
                        nodes[i][j].set_score(topScore)
                        nodes[i][j].set_top(1)
                    else:
                        topScore=int(nodes[i-1][j].get_score())+gapScore
                        nodes[i][j].set_score(topScore)
                        nodes[i][j].set_top(1)
                else:
                    if(isMultipleAlignment==1):
                        elementSeq2 = sequence2[i-1].split(",")
                        leftGapScore = elementSeq2.__len__() * gapScore
                        elementSeq1 = sequence1[j-1].split(",")
                        topGapScore = elementSeq1.__len__() * gapScore

                        leftScore=int(nodes[i][j-1].get_score())+leftGapScore
                        topScore = int(nodes[i - 1][j].get_score()) + topGapScore
                        totalAddingScore=0
                        lenght1=elementSeq1.__len__()
                        lenght2=elementSeq2.__len__()
                        totalAddingScore = 0
                        for a in range(lenght1):
                            for b in range(lenght2):
                                if(elementSeq1[a]=="-"):
                                    elementSeq1[a]="*"
                                if(elementSeq2[b]== "-"):
                                    elementSeq2[b]= "*"
                                totalAddingScore += getScoreFromBlosum62(elementSeq1[a], elementSeq2[b])

                        diagonalScore = int(nodes[i - 1][j-1].get_score())+totalAddingScore
                        maximumScore = max(leftScore, diagonalScore, topScore)
                        nodes[i][j].set_score(maximumScore)
                        if diagonalScore == maximumScore:
                            nodes[i][j].set_diagonal(1)
                        if leftScore == maximumScore:
                            nodes[i][j].set_left(1)
                        if topScore == maximumScore:
                            nodes[i][j].set_top(1)
                    else:
                        leftScore = int(nodes[i][j - 1].get_score()) + gapScore
                        topScore = nodes[i - 1][j].get_score() + gapScore
                        diagonalScore = int(nodes[i - 1][j - 1].get_score()) + getScoreFromBlosum62(sequence1[j - 1],sequence2[i - 1])
                        maximumScore = max(leftScore, diagonalScore, topScore)
                        nodes[i][j].set_score(maximumScore)
                        if diagonalScore == maximumScore:
                            nodes[i][j].set_diagonal(1)
                        if leftScore == maximumScore:
                            nodes[i][j].set_left(1)
                        if topScore == maximumScore:
                            nodes[i][j].set_top(1)
    backtrack(nodes)
    if(isMultipleAlignment==0):
        getMaxAlignmentScore(nodes,sequence1,sequence2)
    else:
        getProfileForMultipleAlignment(nodes,sequence1,sequence2)

def createTempAlignmentTree(): #align edilenlerin hepsini dictionarye ekliyo
    for i in range(0, len(alignmentTree)):
        item = re.split("[(,)]+", alignmentTree[i])
        if (item[1].find("s") > -1 and item[2].find("s") > -1):  # ikisinde de s varsa
            tempAlignmentTree[i] = alignmentTree[i]
        elif (item[1].find("s") <= -1 and item[2].find("s") > -1):  # ilkinde s yok, ikincisinde s varsa
            index = re.split("[()]+", alignmentTree[int(item[1])])
            #tempAlignmentTree[i] = "("+index[1] + "," + item[2]+")"
            tempIndex = re.split("[(,)]+", index[1])
            if(tempIndex[0].find("s") <= -1 and tempIndex[1].find("s") > -1):
                tempIn = re.split("[(,)]+", alignmentTree[int(tempIndex[0])])
                tempAlignmentTree[i] = "(" + tempIn[1]+","+tempIn[2]+"," + tempIndex[1]+ "," + item[2] + ")"
            else:
                tempAlignmentTree[i] = "(" + index[1] + "," + item[2] + ")"
        elif (item[1].find("s") > -1 and item[2].find("s") <= -1):  # ilkinde s var, ikincisinde yoksa
            index = re.split("[()]+", alignmentTree[int(item[2])])
            tempAlignmentTree[i] = "("+item[1] + "," + index[1]+")"
        elif (item[1].find("s") <= -1 and item[2].find("s") <= -1):  # ikisinde de s yoksa
            index1 = re.split("[()]+", alignmentTree[int(item[1])])
            index2 = re.split("[()]+", alignmentTree[int(item[2])])
            tempAlignmentTree[i] ="("+ index1[1] + "," + index2[1]+")"

def updateSimilarityMatrix(firstMatrix):
    global parent
    maxScore = 0
    column = 0
    row = 0
    tempMatrix = [[0 for i in range(len(firstMatrix)-1)] for j in range(len(firstMatrix[0])-1)]
    for i in range(1, len(firstMatrix)):
        for j in range(1, len(firstMatrix[0])):
            if(i < j):
                if(firstMatrix[i][j] >= maxScore):
                    maxScore = firstMatrix[i][j]
                    column = j
                    row = i
    global flag
    flag = 0
    if(parent == 0):
        temp = "(" + firstMatrix[row][0] + "," + firstMatrix[0][column] + ")"
        alignmentTree[parent] = temp
        parent += 1

        item = re.split("[(,)]+", alignmentTree[0])
        if (item[1].find("s") > -1 and item[2].find("s") > -1):  # ikisinde de s varsa
            tempAlignmentTree[0] = alignmentTree[0]
    else:
        temp = ""
        for j in range(0, len(tempAlignmentTree)):
            if((tempAlignmentTree[j].find(firstMatrix[0][column]) > -1 or tempAlignmentTree[j].find(firstMatrix[row][0]) > -1) and flag == 0):
                flag = 1
                temp = "(" + str(j)
            elif(tempAlignmentTree[j].find(firstMatrix[row][0]) > -1 or tempAlignmentTree[j].find(firstMatrix[0][column]) > -1):
                if(flag == 1):
                    temp += ","+str(j) + ")"
                flag = 1
        if(temp.find(")") <= -1 or temp.find("(") <= -1):
            for i in range(0, len(alignmentTree)):
                if (tempAlignmentTree[i].find(firstMatrix[0][column]) > -1):
                    temp = "(" + str(i) + "," + firstMatrix[row][0] + ")"
                elif (tempAlignmentTree[i].find(firstMatrix[row][0]) > -1):
                    temp = "(" + str(i) + "," + firstMatrix[0][column] + ")"
                else:
                    temp = "("+firstMatrix[row][0]+","+firstMatrix[0][column]+")"
                flag = 0
        alignmentTree[parent] = temp
        parent += 1
        createTempAlignmentTree()

    flag = 0
    t = len(firstMatrix)
    i = 0
    j = 0
    while i < t:
        if (i == 0):
            tempMatrix[0][0] = "s0"
            i += 1
            j += 1
        elif (i == row and flag == 0):
            flag = 1
            tempMatrix[0][i] = str(firstMatrix[0][j]) +","+ str(firstMatrix[0][column])
            tempMatrix[i][0] = str(firstMatrix[j][0]) +","+  str(firstMatrix[0][column])
            i += 1
            j += 1
        elif (i != column):
            if (i >= len(tempMatrix) and j < len(firstMatrix)):
                tempMatrix[0][i - 1] = firstMatrix[0][j]
                tempMatrix[i - 1][0] = firstMatrix[j][0]
            elif(j < len(firstMatrix)):
                tempMatrix[0][i] = firstMatrix[0][j]
                tempMatrix[i][0] = firstMatrix[j][0]
            i += 1
            j += 1
        elif(i < len(tempMatrix)):
            if(tempMatrix[i][0] == 0):
                j += 1
                tempMatrix[i][0] = firstMatrix[j][0]
                tempMatrix[0][i] = firstMatrix[0][j]
            i += 1
            j += 1
        else:
            i += 1
            j += 1
    for i in range(1, len(tempMatrix)):
        for j in range(1, len(tempMatrix)):
            if( i >= j ):
                tempMatrix[i][j] = "-"

    tempMatrix = updateScore(column, row, firstMatrix, tempMatrix)
    if(len(tempMatrix) != 2):
        updateSimilarityMatrix(tempMatrix)

def updateScore(maxScoreY, maxScoreX, firstMatrix, tempMatrix):
    x = len(tempMatrix)
    for i in range(1, len(tempMatrix)):
        for j in range(1, len(tempMatrix)):
            if (i < j):
                if(i == maxScoreX):
                    for k in range(1, len(firstMatrix)):
                        if(tempMatrix[0][j] == firstMatrix[k][0]):
                            if(maxScoreY < k):
                                tempMatrix[i][j] = (firstMatrix[maxScoreX][k] + firstMatrix[maxScoreY][k]) / 2
                            else:
                                tempMatrix[i][j] = (firstMatrix[maxScoreX][k] + firstMatrix[k][maxScoreY]) / 2
                else:
                    c = 0
                    b = 0
                    for k in range(1, len(firstMatrix)):
                        if(tempMatrix[i][0] == firstMatrix[k][0]):
                            c = k
                    for k in range(1, len(firstMatrix)):
                        if (tempMatrix[0][j] == firstMatrix[k][0]):
                            b = k
                    if(c != 0 and b != 0):
                        tempMatrix[i][j] = firstMatrix[c][b]
                    else:
                        tempMatrix[i][j] = firstMatrix[i][j]
    return tempMatrix

def getFirstSimilarityMatrix():
    similarityMatrix = [[0 for i in range(sequences.__len__()+1)] for j in range(sequences.__len__()+1)]
    similarityMatrix[0][0] = "s0"
    for i in range(1, sequences.__len__()+1):
        similarityMatrix[0][i] = "s"+str(i)
        similarityMatrix[i][0] = "s" + str(i)
    for i in range(1, sequences.__len__()+1):
        for j in range(1, sequences.__len__() + 1):
            if( i >= j ):
                similarityMatrix[i][j] = "-"
            if (i < j):
                getAlignment(sequences[i-1], sequences[j-1])
                similarityScore = max(maxScore)
                del maxScore[:]
                similarityMatrix[i][j] = round(similarityScore,2)
    updateSimilarityMatrix(similarityMatrix)

getFirstSimilarityMatrix()

def treeCreate(parentNode,key):
    item = re.split("[(,)]+", alignmentTree[int(key)])
    y = Node(item[1], parent=parentNode)
    x = Node(item[2], parent=parentNode)
    if(key == 0):
        return 1
    if (item[1].find("s") <= -1):
        treeCreate(y, item[1])
    if (item[2].find("s") <= -1):
        treeCreate(x, item[2])

item = re.split("[(,)]+", alignmentTree[len(alignmentTree)-1])
rootNode =Node("")
y=Node(item[1], parent = rootNode)
x=Node(item[2], parent = rootNode)
if (item[1].find("s") <= -1):
    treeCreate(y,item[1])
if (item[2].find("s") <= -1):
    treeCreate(x,item[2])

print("GUIDE TREE")
print("------------------")
for pre, fill, node in RenderTree(rootNode):
    if(node.name.isdigit() is True):
        print("%s" % (pre))
    else:
        if(node.name != ""):
            temp = (node.name).split("s")
            for name in range(0, len(sequenceName)):
                if((int(temp[1]) -1) == name):
                    print("%s%s" % (pre, sequenceName[name]))
                    break

multipleAlignment()