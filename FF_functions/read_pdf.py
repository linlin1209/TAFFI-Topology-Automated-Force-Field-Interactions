#!/bin/env python                                                                                                                                                              
# importing required modules 
import sys,PyPDF2 
import argparse,os,datetime,fnmatch,os,re

def main(argv):
  
   # creating a pdf file object 
   pdfFileObj = open('/home/lin1209/pdf/2005_general_drude.pdf', 'rb') 
     
   # creating a pdf reader object 
   pdfReader = PyPDF2.PdfFileReader(pdfFileObj) 
     
   # printing number of pages in pdf file 
   #print(pdfReader.numPages) 
     
   # creating a page object 
   pageObj = pdfReader.getPage(6) 
     
   # extracting text from page 
   text = pageObj.extractText()
   print(text.split('\n'))
   quit()

   with open('scrape.txt','w') as f:
      f.write("{:<40s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s}\n".format("name","formula","exp","ahp Miller",'cc-pVdz','ratio','aug-ccPVdz','ratio')) 
      flag = 0
      for i in text.split('\n'):
         if i == 'waterH': 
            flag = 1
         if flag == 1:
            char=split_str(i)
            for count_j,j in enumerate(char):
               if j.isupper(): 
                  char.pop(count_j)
                  store = j
            name = "".join(char)
            flag = 2 
            continue
         if flag == 2:
            numbers = i.split('.')
            lead = []
            tail = []
            print(numbers)
            for count_j,j in enumerate(numbers):
               tmp = split_str(j)
               if count_j == 0: 
                 lead += tmp[-1] 
                 formula = "".join([store]+tmp[:-1])
                 continue
               if "".join(tmp[2:]) != "":
                  lead.append("".join(tmp[2:]))
               tail.append(tmp[0]+tmp[1])
            print(lead)
            print(tail)
            number = []
            for count_j,j in enumerate(lead):
               number.append(float(".".join([j]+[tail[count_j]])))
            f.write("{:<40s} {:<15s} ".format(name,formula))
            for j in number:
               f.write("{:<15.2f} ".format(j))
            f.write("\n")
            flag =1
            if i == 'AVER0.831.01':break
     
   # closing the pdf file object 
   pdfFileObj.close() 

def split_str(word): 
    return [char for char in word]

if __name__ == "__main__":
   main(sys.argv[1:])
