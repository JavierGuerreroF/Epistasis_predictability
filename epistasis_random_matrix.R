### CLEAN
## RAMDOM MATRIX

#!/usr/bin/env Rscript

## ARGUMENTS
args = commandArgs(trailingOnly=TRUE)

#test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied: seed for runif", call.=FALSE)
} 

t<-args[1]


t<-1

max_g<-c(24,27,30,33,36)
generations<-max(max_g)
N_vec<-c(1,10,100,1000) #Bottleneck
t_ld=23 			# generations in the Luria-Delbruck period
mu=1/(2^t_ld)
replicas<-200

# initial r_vec

r_vec<-c(2.0, 1.2,1.4, 1.6 ,1.8,1.9)
phenotypes<-r_vec[2:length(r_vec)]

# kvalues for the formula
cou=0


clock <- Sys.time()
print('start time:')
print(clock)


set.seed(t)
epistasis_matrix<-matrix(runif(25,0,2.2),5,5)
rownames(epistasis_matrix)<-phenotypes
colnames(epistasis_matrix)<-phenotypes
cat("runif seed = ", t,"\n")

fig <- plot_ly(z = epistasis_matrix,
               x = phenotypes,
               y = phenotypes)
fig <- fig %>% add_surface()

fig


####
x = c(1,2,3,4,5)
y = c(1,2,3,4,5)
z = epistasis_matrix


fig <- plot_ly(
  type = 'surface',
  contours = list(
    z = list(show = TRUE, start = 0, end = 2.2, size = 0.1)),
  x = ~x,
  y = ~y,
  z = ~z,)
fig <- fig %>% layout(
  scene = list(
    xaxis = list(nticks = 20),
    zaxis = list(nticks = 4)))
    #camera = list(eye = list(x = 0, y = -1, z = 0.5)),
    #aspectratio = list(x = .9, y = .8, z = 0.2)))

fig

###
for (N in N_vec){
  
  print(clock)
  print(epistasis_matrix)
  
  cat('Bottleneck = ',N,'\n')

##### POP CREATION

x<-matrix(0,max(max_g),length(r_vec))




r_vec_mod<-0
sample_N<-0

for (gen1 in max_g){
  for (gen2 in max_g){
    
    if(gen1==gen2){
      
   
    
    final_table<-matrix(0,5,5)
    cat('\n',gen1," ",gen2,"\n")
    
    for (repl in 1:replicas){
      
      #print('Replica')
      #print(n)
      
      if (repl%%(replicas/10)==0){
        cat(' ',repl)
      }
      
      #cat(i," ",j,"\n")
      gen<-c(gen1,gen2)
      
      result<-matrix(0,N,2)
      
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
              c=c+1
              
              
              # .. reinitialization
              x[]=0
              
              mut_tot=0 
              x[1,1]=1	
              
              
              for (g in 1:(generations-1)){
                #cat(g,"\n")
                #cat(".")
                
                
                
                
                
                #.. population growth
                
                
                #epistasis growth
                #if (n==1){
                x[g+1,]=r_vec*x[g,]
                #  r_vec_mod<-0
                #}else{
                #  r_vec_mod<-append(phenotypes[sample_N],epistasis_matrix[sample_N,])
                #  x[g+1,]=r_vec_mod*x[g,]
                #}
                
                #x[g+1,]=r_vec*x[g,] # growth
                
                
                
                
                
                mut <- rpois(1,lambda=(x[g+1,1]-x[g+1,1]/2)*mu)  	# mutations per division!
                #cat('mutations =',mut,'\n')
                x[g+1,1]=x[g+1,1] - mut				# update of population matrix
                mut_tot=mut_tot+mut
                
                
                x[x<=0]<-0					
                if (mut>0) { 					# allocation of mutations
                  for (z in 1:mut) {
                    
                    clone<-sample(seq(2,6,1),1)
                    
                    
                    
                    x[g+1,clone]=x[g+1,clone]+1
                  }  	
                  
                }
                
                
                
                
              }
              #print(r_vec_mod)
              
              if (sum(x[w,2:length(r_vec)]>0)){
                pick<-1
                #print('No mutation found')
                #print(c)
                
              }
              
              
              ## crear aqui vector/matrix prob
              
              
              #cat("\n")
            }
            
            
            prob<-x[w,2:length(r_vec)]/sum(x[w,2:length(r_vec)])
            
            
            #result[n]<-sample(r_vec[2:length(r_vec)],1,prob=prob)
            #result[n]<-sample(seq(1,5,1),1,prob=prob)
            
            #print(r_vec_mod)
            
            sample_N<-sample(seq(1,5,1),N,prob=prob,replace=TRUE)
            #print(sample_N)
            #result[,n]<-sample_N
            
            
          }
          
          ###### POP GEN 1 #############################################
          
          x_first<-x
        }else{
          
          ###### POP GEN 2 #############################################
          table_data<-as.data.frame(table(sample_N))
          sample_data<-matrix(0,5,2)
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
              
              
              
              #print(ph)
              #print(phenotypes[ph])
              ##################### POP gen
              mut_tot=0
              while (mut_tot<1) {
                
                pick=0
                c=0
                while (pick<1){
                  c=c+1
                  
                  
                  # .. reinitialization
                  x[]=0
                  
                  mut_tot=0 
                  ## modificar la x innical de acuerdo a la tabla de frecuencias del sample_N  
                  
                  x[1,1]=sample_data[ph,2]	
                  
                  
                  for (g in 1:(generations-1)){
                    #cat(g,"\n")
                    #cat(".")
                    
                    
                    
                    
                    
                    #.. population growth
                    
                    
                    #epistasis growth
                    
                    r_vec_mod<-0
                    r_vec_mod<-append(phenotypes[ph],epistasis_matrix[ph,])
                    x[g+1,]=r_vec_mod*x[g,]
                    
                    
                    #x[g+1,]=r_vec*x[g,] # growth
                    
                    
                    
                    
                    
                    mut <- rpois(1,lambda=(x[g+1,1]-x[g+1,1]/2)*mu)  	# mutations per division!
                    #cat('mutations =',mut,'\n')
                    x[g+1,1]=x[g+1,1] - mut				# update of population matrix
                    mut_tot=mut_tot+mut
                    
                    
                    x[x<=0]<-0					
                    if (mut>0) { 					# allocation of mutations
                      for (z in 1:mut) {
                        
                        clone<-sample(seq(2,6,1),1)
                        
                        
                        
                        x[g+1,clone]=x[g+1,clone]+1
                      }  	
                      
                    }
                    
                    
                    
                    
                  }
                  #print(r_vec_mod)
                  
                  if (sum(x[w,2:length(r_vec)]>0)){
                    pick<-1
                    #print('No mutation found')
                    #print('POP2')
                    #print(c)
                    
                  }
                  
                  
                  ## crear aqui vector/matrix prob
                  
                  
                  #cat("\n")
                }
                
                
                #prob<-x[w,2:length(r_vec)]/sum(x[w,2:length(r_vec)])
                
                
                #result[n]<-sample(r_vec[2:length(r_vec)],1,prob=prob)
                #result[n]<-sample(seq(1,5,1),1,prob=prob)
                
                #print(r_vec_mod)
                
                #sample_N_final<-sample(seq(1,phenotypes**2,1),N,prob=prob,replace=TRUE)
                #print(sample_N)
                #result[,n]<-sample_N
                
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
          
          sample_N_final<-sample(seq(1,length(phenotypes)**2,1),N,prob=prob,replace=TRUE)
          sample_N_final_vec<-rep(0,25)
          for(sample in sample_N_final){
            sample_N_final_vec[sample]<-sample_N_final_vec[sample]+1
          }
          
          ###### POP GEN 2 #############################################
          
          
          
        }
        
        
        
      }
      
      #print(result)
      result_matrix<-matrix(sample_N_final_vec,length(phenotypes),length(phenotypes),byrow=TRUE)
      
      final_table<-final_table+result_matrix
      
      #final_table[result[1],result[2]]<-final_table[result[1],result[2]]+1
      
      
    }
    cat("\n")
    name<-paste0("data_",gen1,"_",gen2,"_rep_",replicas,"_bn_",N,"_seed_",t,".dat")
    colnames(final_table)<-r_vec[2:length(r_vec)]
    rownames(final_table)<-r_vec[2:length(r_vec)]
    write.table(final_table, file = name, sep = " ")
    #print(sum(final_table))
    print(name)
    print(Sys.time())
    
  }
    
    
  }
  
  
  
  
}


}
  

print(Sys.time()-clock)

## END OF SIMULATION

