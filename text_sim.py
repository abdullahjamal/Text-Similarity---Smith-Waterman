'''
To compile code, run python prog1_new.py --src=pepper-src.txt --target=pepper-tgt.txt
You just need to pass the source text file to src and target text file to target
'''
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Text Similarity inputs')
parser.add_argument('--src', type=str,
                   help='source file input')
parser.add_argument('--target', type=str,
                   help='input target file')

args = parser.parse_args()

gap_insert = -1
gap_delete = -1
mismatch = -1
match = 2

print "\nText Similarity Analysis\n"

print"Source file:"+' '+args.src
print"Target file:"+' '+args.target


''' 
Below function takes the file as an argument 
and convert the letters in to lower case and split
the string w.r.t to space. It then applies the second and third rule.
'''
def normalize_tokens(input_file):

    tokens = input_file.lower()
    tokens = tokens.split()
    num_tokens = []

    extra_tokens = []
    check = False
    rule2 = False
    rule3 = False

    # checking if token is alphanumeric
    for i in range(len(tokens)):
        token = tokens[i]
        for j in range(len(token)):
            if token[j].isalnum():
                check = True

        # check if first char is non-alphanumeric
        ptr = 0
        if check == True and token[0].isalnum() == False:
            rule2 = True

            while (token[ptr].isalnum() == False):
                extra_tokens.append(token[ptr])
                ptr += 1
            slashed_token = token[ptr:len(token)]

            for k in range(len(extra_tokens)):
                num_tokens.append(extra_tokens[k])

            if token[-1].isalnum() == True:
                num_tokens.append(slashed_token)
  
        extra_tokens = []
        if check == True and token[-1].isalnum() == False:
            rule3 = True
            if rule2 == True:
                token = token[2-1:len(token)]
            ptr = -1
            while (token[ptr].isalnum() == False):
                extra_tokens.append(token[ptr])
                ptr -= 1
            slashed_token = token[0:len(token)+ptr+1]
            num_tokens.append(slashed_token)

            for k in range (len(extra_tokens)):
                num_tokens.append(extra_tokens[len(extra_tokens)-1-k])

        else:
            if check == False or (rule2 == False and rule3 == False):
                num_tokens.append(token)

        check = False
        rule2 = False
        rule3 = False
        extra_tokens = []
        ptr = 0

    return num_tokens

''' 
Below function takes the list of tokens and apply the 
last rule. It checks the condition for apostrophe (''s'), 'not'
,and 'am' 
'''
def final_norm_lastrule(num_tokens):
        
    aps_list = []
    for ind,wd in enumerate(num_tokens):
        if wd[-2:] == "'s":
            aps_list.append(wd[:-2])
            aps_list.append("'s")
        else:
            aps_list.append(wd)
        
    not_list = []
    for ind,wd in enumerate(aps_list):
        
        if wd[-3:] == "n't":
           not_list.append(wd[:-3])
           not_list.append("not")
        else:
           not_list.append(wd)
        
    am_list = []
    for ind,wd in enumerate(not_list):
        
        if wd[-2:] == "'m":
           am_list.append(wd[:-2])
           am_list.append("am")
        else:
           am_list.append(wd)
        
    
    return am_list
    
    
f = open(args.src,'r')
src_str = f.read()
f.close()

ft = open(args.target,'r')
tgt_str = ft.read()
ft.close()

src_norm = normalize_tokens(src_str)
tgt_norm = normalize_tokens(tgt_str)

src_norm_new = final_norm_lastrule(src_norm)
tgt_norm_new = final_norm_lastrule(tgt_norm)

print"\nRaw Tokens:\n"
print'     '+"Source >"+' '+src_str.lower()
print'     '+"Target >"+' '+tgt_str.lower()

print"\nNormalized Tokens:\n"
print'     '+"Source >"+' '+' '.join(src_norm_new)
print'     '+"Target >"+' '+' '.join(tgt_norm_new)

src_norm_new.insert(0,'#')
tgt_norm_new.insert(0,'#')

rows = len(src_norm_new) 
cols = len(tgt_norm_new) 

matrix = np.zeros((rows,cols))
backtable = np.zeros((rows,cols))
backtable_cp = np.chararray((rows,cols),itemsize=2)
backtable_cp[0,:] = ' '
backtable_cp[:,0] = ' '


''' 
Function below generate the two strings as source and 
target by using backtrace table
'''
def trace(backtable,i,j):
    align1 = ""
    align2 = ""
    while backtable[i][j] != 0:
      if backtable[i][j] == 3:
        align1 = align1+" "+src_norm_new[i]
        align2 = align2+" "+tgt_norm_new[j]
        i-=1
        j-=1
      elif backtable[i][j] == 2:
        align1 = align1+" "+'-'
        align2 = align2+" "+tgt_norm_new[j]
        j-=1
      elif backtable[i][j] == 1:
        align1 = align1+" "+src_norm_new[i]
        align2 = align2+" "+'-'
        i-=1
    
    return ' '.join(reversed(align1.split(' '))),' '.join(reversed(align2.split(' '))),i,j
    
    
'''
The function below will check both strings and
assign symbols like 'i','d','-'
'''
def alignment_string(str1,str2):

    align_symbol = []
    for b1, b2 in zip(str1.split(' '), str2.split(' ')):
        if b1 == b2:
            align_symbol.append(' ')
            
        elif '-' in b2:
            align_symbol.append('d')
            
        elif '-' in b1:
            align_symbol.append('i')
            
        else:
            align_symbol.append('s')
           

    return align_symbol #' '.join(align_symbol) if I uncomment the line, then it outputs the symbols in correct format for only gene file. That's why the function 
                                                # returns only list of symbols. The format is not corrrect, but the symbols are correct for each test case.
    
    
''' 
The function calculates score for each index for
edit distance table
'''
def score_matrix(matrix, x, y):

    similarity = match if src_norm_new[x] == tgt_norm_new[y] else mismatch

    diag_score = matrix[x - 1][y - 1] + similarity
    up_score   = matrix[x - 1][y] + gap_insert
    left_score = matrix[x][y - 1] + gap_insert

    return max(0, diag_score, up_score, left_score),diag_score,up_score,left_score
    
''' calculate edit distance table '''
max_score = 0
max_pos  = None
for i in range(1,rows):
    for j in range(1,cols):
        score,dg_scr,up_scr,lt_scr = score_matrix(matrix, i, j)
        if score == 0:
           backtable[i][j] = 0
           backtable_cp[i][j] = ' '
        elif score == lt_scr:
            backtable[i][j] = 2
            backtable_cp[i][j] ='LT'
        elif score == up_scr:
            backtable[i][j] = 1
            backtable_cp[i][j] = 'UP'
        elif score == dg_scr:
            backtable[i][j] = 3
            backtable_cp[i][j] = 'DI'
        if score > max_score:
            max_score = score
            max_pos = (i,j)
        matrix[i][j] = score


''' creating index for both edit distance and backtrace table'''
new_src = [src_norm_new[i][0:3] for i in range(0,len(src_norm_new))]
new_tgt = [tgt_norm_new[i][0:3] for i in range(0,len(tgt_norm_new))]

src_list = [s for s in range(0,len(new_src))]
tgt_list = [t for t in range(0,len(new_tgt))]

row_list =[]
row_list.append(src_list)
row_list.append(new_src)

col_list =[]
col_list.append(tgt_list)
col_list.append(new_tgt)

''' matrix array contains edit distance table'''

print "\nEdit Distance Table:\n"
df_dist = pd.DataFrame(data=matrix,
                       index=row_list,
                       columns=col_list)

print df_dist

''' Backtrace Table - Alignment'''
''' backtable_cp array contains backtrace table'''

print "\nBacktrace Table:\n"

df_bt = pd.DataFrame(data=backtable_cp,
                     index=row_list,
                     columns=col_list)

print df_bt


max_val = np.amax(matrix)

print "\nMaximum value in distance table:",int(max_val)


ind_max = np.where(matrix==max_val)
ind1_max = zip(ind_max[0],ind_max[1])

print "\nMaxima:"
for m in range(0,len(ind1_max)):
    print "    ",list(ind1_max[m]),"\n"

print"\nMaximal-similarity alignments:\n"


for k in range(0,len(ind1_max)):
    i,j = ind1_max[k][0],ind1_max[k][1]
    
    str1,str2,a,b = trace(backtable,i,j)
    print "    Alignment"+" "+str(k)+" (length"+"",len(str1.split()),"):"
    print "       Source at   "+str(a)+": ",str1
    print "       Target at   "+str(b)+": ",str2

    symbols = alignment_string(str1,str2)
    
    print "       Edit action  : ",symbols
