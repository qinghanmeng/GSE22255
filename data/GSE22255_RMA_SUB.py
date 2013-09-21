#generate a dictionary.

f_dic = open("/home/qinghan/workspace/GSE22255/GPL/GPL570-13270-ID-TITLE-SYMBOL")
dic = []
for line in f_dic:
	words = line.strip("\n").split("\t")
	dic.append(words[0])


f = open("GSE22255_RMA.TXT")
f_sub = open("GSE22255_RMA_SUB.TXT","w") 
for line in f:
	words = line.split("\t")
	if (words[0] in dic):
		f_sub.write(line)


