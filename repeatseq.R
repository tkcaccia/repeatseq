# repeatqeq is a function that screens DNA sequences for interspersed repeats 
# and low complexity DNA sequences. The output of the program is a detailed annotation 
# of the repeats that are present in the query sequences as well as the number of
# duplicates in the query sequences.


repeatseq = function(str,len=c(1,2,3,4,5)){
  str=as.character(str)
  lista=strsplit(str,"")
  tot=length(lista)
  n_len=length(len)
  dupl_matrix=matrix(0,nrow=tot,ncol=n_len)

  for(h in 1:n_len){
    for(j in 1:tot){
      gene=lista[[j]]
      n=length(gene)-len[h]+1
      dupl=NULL
      for(i in 1:n){
        a=i
        a1=a+len[h]-1
        primer=gene[a:a1]
        dupl[i]=0
        while(a<=n & a>0){
          if(all(gene[a:a1]==primer)){
            dupl[i]=dupl[i]+1
            a=a+len[h]
            a1=a1+len[h]
          }else{
            a=-1
          }
        }
      }
      dupl_matrix[j,h]=max(dupl)
    }
  }
  colnames(dupl_matrix)=as.vector(len)
  rownames(dupl_matrix)=str
  t=table(str)
  duplicated=t[str]
  dupl_matrix=cbind(dupl_matrix,duplicated)
  dupl_matrix
}


str=c("CGCGCGCGCGCGCGCGCGCGCGATGTACTT",
      "TGACCTCTGGAATCGAGGTCAG",
      "TGTTGTTGTTGTTGTTGTTGTTGTTGTTGT",
      "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
      "CCTTCCTTCCTTCCTTCCTTCCTTCCTT",
      "ACGTTGACGGACCTAACTGACTGCCGTGACCTTAC",
      "ACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTA",
      "AAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTT",
      "TGAACCTATTAGGCCAATTTTAACCCGTACCAATGCCGGGGT",
      "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
      "AAGTGCGCGGATCAC",
      "AAGTGCGCGGATCAC")

repeatseq(str)
