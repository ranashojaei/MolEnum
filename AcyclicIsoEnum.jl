

using DataFrames 
using Plots 
using LaTeXStrings #to use latex in my plot labaling


function count_acyclic_isomers(L:: Integer, M:: Integer, N:: Integer) 


    b=zeros(BigInt,L,M,N) 
    a=zeros(BigInt,L,M,N) 

    ############## Initial Conditions ###############
    a[1,1,1]=0;
    a[2,1,1]=0; 
    a[1,2,1]=0;  
    a[1,1,2]= a[1,2,2]=1; #term A0 and A1 in the figure in the paper (h+oh)
    


    b[1,1,1]=0; 
    b[2,1,1]=0;
    b[1,2,1]=1;  #term B1 in the figure in the paper (o)
    b[1,1,2]=0; 


    results = []   # for saving the number of isomers to csv file
    ################################################
    for l=2:L 
        for m=1:M 
            for n=1:N 
                ################### coefficients b ################## 

                ### Term B2 
                TermB2=0 
                S=l-2 
                Sp=m-1 
                Spp=n-1 

                for s=0:S 
                    for sp=0:Sp 
                        for spp=0:Spp 
                            TermB2 += a[s+1,sp+1,spp+1] * a[l-s-1,m-sp,n-spp]
                        end 
                    end 
                end 

            
                b[l,m,n]= TermB2    


                ################### coefficients a ################## 
                #####################################################

                ### Term A3 
                TermA3=0 
                S= l-2 
                Sp= m-1 
                Spp= n-1 

                for s=0:S 
                    for sp=0:Sp 
                        for spp=0:Spp 

                            TermA3 += a[l-s-1,m-sp,n-spp] * b[s+1,sp+1,spp+1] 
                        end
                    end
                end 
                

                ### Term A2_1 
                TermA2_1=0 
                S= l-2 
                Sp= m-1 
                Spp= n-1 
                

                for s=0:S 
                    for sp=0:Sp 
                        for spp=0:Spp 

                            T=l-s-2 
                            Tp=m-sp-1
                            Tpp=n-spp-1 

                            for t=0:T 
                                for tp=0:Tp 
                                    for tpp=0:Tpp 

                                        TermA2_1 += a[t+1,tp+1,tpp+1] * a[s+1,sp+1,spp+1] * a[l-s-t-1,m-sp-tp,n-spp-tpp] 

                                    end 
                                end 
                            end 

                        end 
                    end 
                end 


                ### Term A2_2
                TermA2_2=0 
                if (l+1)%3==0 && (m+2)%3==0 && (n+2)%3==0 
                    TermA2_2= a[(l+1)÷3, (m+2)÷3, (n+2)÷3]
                end 
                
    
                ### Term A4
                TermA4=0 
                if l>=3
                TermA4 = a[l-2, m,n] 
                end 


                # the number of rooted acyclic isomers
                a[l,m,n]= TermA4 + (2*TermA2_2 + TermA2_1)÷3 + TermA3   
                

                #to store the number in a dataframe/csv file
                value = a[l, m, n]  
                if value != 0
                    push!(results, (l-1, m-1, n-1, value))  # store (l, m, n, isomer)
                    println(l-1, ",", m-1, ",", n-1, " : ", value)
                end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

            end 
        end 
    end 
    Co_a = DataFrame(results, [:l, :m, :n, :isomer])  #Coefficients a
    return Co_a
end  




#Enter L carbon, M oxygen, and N Hydrogen atom in count_acyclic_isomers(L,M,N) 
Co_a = count_acyclic_isomers(10,6,21)  #typeof(Co_a): DataFrame
Iso_List=Matrix(Co_a)  #Matrix



# Fig.6: Visualization of isomer count vs. C, O, and H atoms
# For reproducing Fig.6 in the paper: set L = 21, M = 16, N = 44 , to ensure sufficient isomer generation
#########################################################
# Fig 6. (a) in the paper for O=10, H=5 
carbonNumber=[]
isoNumber_C=[] 
for i=1:size(Iso_List, 1)
    if Iso_List[i,2]==10 && Iso_List[i,3]==5   
        push!(carbonNumber, (Iso_List[i,1]))
        push!(isoNumber_C, (Iso_List[i,4]))
    end
end 
carbonNumber
isoNumber_C


Plot1= scatter(
    log.(carbonNumber), log.(isoNumber_C) ,markersize=6,   
    grid=false, markercolor=:steelblue, 
    size=(600,600) ,legend=false, 
    xlabel=L"\ln(\mathrm{Number\ of\ Carbon\ Atoms})",  
    ylabel=L"\ln(\mathrm{Number\ of\ Isomers})",  
    framestyle = :box 
)   


#########################################################
# Fig 6. (b) in the paper for C=10, H=15 (must be odd) 
OxygenNumber=[]
isoNumber_O=[] 
for i=1:size(Iso_List, 1)
    if Iso_List[i,1]==10 && Iso_List[i,3]==15 && Iso_List[i,2]!=0   #Iso_List[i,3] must be odd, up to N

        push!(OxygenNumber, (Iso_List[i,2]))
        push!(isoNumber_O, (Iso_List[i,4]))
        
    end
end 
OxygenNumber
isoNumber_O


Plot2=plot(
    OxygenNumber, isoNumber_O, 
    seriestype = :scatter, size = (600, 600),
    markersize=6, legend = false, grid = false,
    xlabel=L"\mathrm{Number\ of\ Oxygen\ Atoms}",  
    ylabel=L"\mathrm{Number\ of\ Isomers}", 
    markercolor=:steelblue, framestyle = :box 
)


#######################################################
# Fig 6. (c) in the paper for C=20, O=15 
HydNum=[]   
isoNumber_H=[] 
for i=1:size(Iso_List, 1)
    if Iso_List[i,1]==20  && Iso_List[i,2]==15
        push!(HydNum, (Iso_List[i,3]))
        push!(isoNumber_H, (Iso_List[i,4]))
        
    end
end 
HydNum
isoNumber_H


Plot3 = plot(
    HydNum, isoNumber_H,
    seriestype = :scatter, size = (600, 600),
    markersize=6, legend = false, grid = false,
    xlabel = L"\mathrm{Number\ of\ Hydrogen\ Atoms}",
    ylabel = L"\mathrm{Number\ of\ Isomers}",
    color=:steelblue,
    framestyle = :box
) 







