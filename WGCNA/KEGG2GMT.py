#!/usr/bin/env python
# coding: utf-8

# In[46]:

import sys
#import pandas

#intable = sys.argv[1]
#outable=sys.argv[2]
# In[115]:


with open('KEGG_annotation/LbrM2903.keg', 'r') as infile:
    data = infile.readlines()


# In[116]:


data = data[3:-4]


# In[117]:


funB = []
for i in data:
    if 'B  ' in i:
#         print(i)
        funB.append(i[9:-1])


# In[118]:


new_data = []
for i in data:
    if 'D  ' in i or 'C  ' in i:
#         print(i)
        new_data.append(i)


# In[119]:


indexes = []
for index, i in enumerate(new_data):
    if 'C    ' in i:
        indexes.append(index)


# In[120]:


indexes


# In[121]:


new_data[0]


# In[122]:


new_data[36]


# In[142]:


tmp = ' '.join([ x.split(';')[0].split(' ')[-1] for x in new_data[1:36]])
tmp


# In[124]:


new = ' '.join(new_data[0].split('  ')[-1].split(' ')[1:])[:-1] + '\t' + tmp
new


# In[125]:


with open('LbraM2903.gmt', 'w') as outfile:
    for i in new:
        outfile.write(i)


# In[154]:


parsed = []
for i in range(len(indexes)-1):
#     print(indexes[i]+1, indexes[i+1]-1)
    parsed.append(' '.join(new_data[indexes[i]].split('  ')[-1].split(' ')[1:])[:-1] + '\tstring\t' + 
                 '\t'.join([ x.split(';')[0].split(' ')[-1] for x in new_data[indexes[i]+1:indexes[i+1]]]) + '\n')
    
with open('LbraM2903.gmt', 'w') as outfile:
    for i in parsed:
        outfile.write(i)


# In[ ]:



