
## 3 CASES MATRICES
# Code used to simulated the cases based on the mathematical model

#!/usr/bin/env Rscript


## ARGUMENTS
args = commandArgs(trailingOnly=TRUE)

#test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied: seed for runif", call.=FALSE)
} 

t<-args[1]

# Parameters

max_g<-c(24,27,30,33,36) # Generations 
generations<-max(max_g)
N_vec<-c(1,10,100,1000)  # Bottleneck
t_ld=23 			           # Generations in the Luria-Delbruck period
mu=1/(2^t_ld)            # Mutation rate
replicas<-200

## SLURM PARAMETERS
param_matrix<-matrix(0,24,2)
param_matrix[,1]<-c(1,1,1,1,1,1,1,1,
                    2,2,2,2,2,2,2,2,
                    3,3,3,3,3,3,3,3)
param_matrix[,2]<-c(-0.6,-0.3,0,0.4,0.8,1,1.3,1.5,
                    -0.6,-0.3,0,0.4,0.8,1,1.3,1.5,
                    -0.6,-0.3,0,0.4,0.8,1,1.3,1.5)

colnames(param_matrix)<-c('case','k2')
rownames(param_matrix)<-seq(1,48,1)

r_vec_list<-list(c(2.0, 1.2,1.3, 1.4 ,1.5,1.6),
                 c(2.0, 1.6,1.7,1.8,1.9,2),
                 c(2.0, 1.2,1.4,1.6,1.8,2))

# case
case<-param_matrix[t,1]
# initial r_vec
r_vec<-r_vec_list[[case]]
phenotypes<-r_vec[2:length(r_vec)]
# kvalues for the formula
k1<-1
k2<-param_matrix[t,2]


## EPISTASIS MATRIX GENERATION
epistasis_matrix<-matrix(0,length(r_vec)-1,length(r_vec)-1)
colnames(epistasis_matrix)<-r_vec[2:length(r_vec)]
rownames(epistasis_matrix)<-r_vec[2:length(r_vec)]

for (i in 1:nrow(epistasis_matrix)){
  for (j in 1:ncol(epistasis_matrix)){
    c1<-2-(phenotypes[i])
    c2<-2-(phenotypes[j])
    
    epistasis_matrix[i,j]<-2-(c2*k2+c1*k1)
  }
}



epistasis_matrix<-round(epistasis_matrix,2)

print(epistasis_matrix)


clock <- Sys.time()
print('start time:')
print(clock)

for (N in N_vec){
  
  print(clock)
  cat('Bottleneck ',N,'\n')

##### POP CREATION

x<-matrix(0,max(max_g),length(r_vec))




r_vec_mod<-0
sample_N<-0

for (gen1 in max_g){
 
    gen2<-gen1
        
    final_table<-matrix(0,nrow(epistasis_matrix),ncol(epistasis_matrix))
    cat('\nGenerations 1st = ',gen1," 2nd = ",gen2,"\n")
    
    for (repl in 1:replicas){
      
   
      if (repl%%(replicas/10)==0){
        cat('.')
      }
      cat(repl," ")
 
      gen<-c(gen1,gen2)
      

      
      n=0
      
      for (w in gen){
        
        n=n+1
        if (n==1){
          ###### POP GEN 1 #############################################
          mut_tot=0
          while (mut_tot<1) {
            
            pick=0
            c=0
            while (pick<1){
              
              # do this till we pick at least one bacteria
              c=c+1
              
              
              # reinitialization
              x[]=0
              
              mut_tot=0 
              x[1,1]=1	
              
              
              for (g in 1:(generations-1)){
               
                # epistasis growth
             
                x[g+1,]=r_vec*x[g,]
             
                
                
                mut <- rpois(1,lambda=(x[g+1,1]-x[g+1,1]/2)*mu)  	# mutations per division
              
                x[g+1,1]=x[g+1,1] - mut				# update of population matrix
                mut_tot=mut_tot+mut
                
                
                x[x<=0]<-0
              
                if (mut>0) { 	
                  
                  # allocation of mutations
                  count<-table(sample(seq(2,length(r_vec),1),mut,replace=TRUE))
                  names<-dimnames(count)[[1]]
                  for (name in names){
                    num<-as.integer(name)
                    x[g+1,num]<-x[g+1,num]+count[name]
                  }
                  
                }
                
                
                
                
              }
           
              if (sum(x[w,2:length(r_vec)]>0)){
                pick<-1
                ## This is to check that there is at least a mutant bacteria in 
                # in the population for the generation selected.
              }
              
              
      
          
            }
            
            
            prob<-x[w,2:length(r_vec)]/sum(x[w,2:length(r_vec)])
            
            
            sample_N<-sample(seq(1,length(phenotypes),1),N,prob=prob,replace=TRUE)
           
            
            
          }
          
          ###### POP GEN 1 #############################################
          
          x_first<-x
        }else{
          
          ###### POP GEN 2 #############################################
          table_data<-as.data.frame(table(sample_N))
          sample_data<-matrix(0,length(phenotypes),2)
          sample_data[,1]<-seq(1,length(phenotypes),1)
         
          for (row in 1:length(phenotypes)){
            
            for (o in 1:nrow(table_data)){
              
              if (row==table_data[o,1]){
                sample_data[row,2]<-table_data[o,2]
              }
            }
          }
          
          
          
          for (ph in 1:length(phenotypes)){
            if (ph ==1){
              final_x_total<-x_first
            }
            if (sample_data[ph,2]>0){
              
              ##################### POP gen
              mut_tot=0
              while (mut_tot<1) {
                
                pick=0
                c=0
                while (pick<1){
                  c=c+1
                  
                  
                  # reinitialization
                  x[]=0
                  
                  mut_tot=0 
                 
                  x[1,1]=sample_data[ph,2]	
                  
                  
                  for (g in 1:(generations-1)){

                    # epistasis growth
                    
                    r_vec_mod<-0
                    r_vec_mod<-append(phenotypes[ph],epistasis_matrix[ph,])
                    x[g+1,]=r_vec_mod*x[g,]
                    
               
     
                    
                    mut <- rpois(1,lambda=(x[g+1,1]-x[g+1,1]/2)*mu)  	# mutations per division!
               
                    x[g+1,1]=x[g+1,1] - mut				# update of population matrix
                    mut_tot=mut_tot+mut
                    
                    
                    x[x<=0]<-0					
                    if (mut>0) { 	
                      
                      # allocation of mutations
                      count<-table(sample(seq(2,length(r_vec),1),mut,replace=TRUE))
                      names<-dimnames(count)[[1]]
                      for (name in names){
                        num<-as.integer(name)
                        x[g+1,num]<-x[g+1,num]+count[name]
                      }
                      
                    }
                    
                    
                    
                    
                    
                  }
                  
        
                  if (sum(x[w,2:length(r_vec)]>0)){
                    pick<-1
                 
                    
                  }
                  
              
                }
                
                
           
                final_x_total<-cbind(final_x_total,x)
              }
            }else{
              

              x[]=0
              final_x_total<-cbind(final_x_total,x)}
            
            #############################
 
          }
          final_dblmut_pop<-c()
          
          for (i in (length(r_vec)+1):ncol(final_x_total)){
            if ((i-1)%%length(r_vec)!=0){
              final_dblmut_pop<-append(final_dblmut_pop,final_x_total[w,i])
              
            }
          }
          prob<-final_dblmut_pop/sum(final_dblmut_pop)
          
          sample_N_final<-sample(seq(1,length(prob),1),N,prob=prob,replace=TRUE)
          sample_N_final_vec<-rep(0,length(prob))
          
          count<-table(sample_N_final)
          names<-dimnames(count)[[1]]
          for (name in names){
            num<-as.integer(name)
            sample_N_final_vec[num]<-sample_N_final_vec[num]+count[name]
          }
          
          
          ###### POP GEN 2 #############################################
          
          
          
        }
        
        
        
      }
      
     
      result_matrix<-matrix(sample_N_final_vec,length(phenotypes),length(phenotypes),byrow=TRUE)
      
      final_table<-final_table+result_matrix

 
    }
    cat("\n")
    ## DATA NAMING
    name<-paste0("data_",gen1,"_",gen2,"_rep_",replicas,"_bn_",
                 N,"_k2value_",k2,"_CASE_",case,".dat")
    colnames(final_table)<-r_vec[2:length(r_vec)]
    rownames(final_table)<-r_vec[2:length(r_vec)]
    write.table(final_table, file = name, sep = " ")
  
    print(name)
    print(Sys.time())
  
}


}
  

print(Sys.time()-clock)

## END OF SIMULATION

