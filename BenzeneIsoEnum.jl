

################# first calculation of coeffs. a #################

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

                if a[l, m, n] != 0
                    println(l-1, ",", m-1, ",", n-1, " : ", a[l,m,n])
                end
                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            

            end 
        end 
    end 

    #Coefficients a
    return a
end  

#Enter L carbon, M oxygen, and N Hydrogen atom in count_acyclic_isomers(L,M,N)    #(  N_max <  ( (2*L)+2 )   )
a = count_acyclic_isomers(11,6,24); #typeof(a): Array{BigInt, 3} 


####################### Benzene ############################
####################### Benzene ############################ 
####################### Benzene ############################

#Enter L carbon, M oxygen, and N Hydrogen atom
L=4     #L_min = 1
M=2     #M_min = 1
N= 9    #N_min = 7

d=zeros(BigInt,L,M,N);    #coeffiecients to reperesent the number of Benzene


for l=1:L 
    for m=1:M 
        for n=1:N 


            ############ One ###############
            one=0

            for s=0:l-1 
                for sp=0:m-1
                    for spp=0:n-1  
                        
                        for t=0:l-s-1 
                            for tp=0:m-sp-1 
                                for tpp=0:n-spp-1 

                                    for u=0:l-s-t-1
                                        for up=0:m-sp-tp-1 
                                            for upp=0:n-spp-tpp-1 

                                                for v=0:l-s-t-u-1
                                                    for vp=0:m-sp-tp-up-1 
                                                        for vpp=0:n-spp-tpp-upp-1 

                                                            for w=0:l-s-t-u-v-1
                                                                for wp=0: m-sp-tp-up-vp-1 
                                                                    for wpp=0:n-spp-tpp-upp-vpp-1 

                                                                        one +=  a[s+1,sp+1,spp+1] * a[t+1,tp+1,tpp+1] * a[u+1,up+1,upp+1] * a[v+1,vp+1,vpp+1] * 
                                                                        a[w+1,wp+1,wpp+1] * a[l-s-t-u-v-w,m-sp-tp-up-vp-wp,n-spp-tpp-upp-vpp-wpp]  

                                                                    end 
                                                                end 
                                                            end 



                                                        end 
                                                    end 
                                                end  


                                            end 
                                        end 
                                    end 


                                end 
                            end 
                        end 



                    end 
                end 
            end 


            ######################### two ############
            two=0 
            if (l+5)%6==0 && (m+5)%6==0 && (n+5)%6==0 
                two= a[(l+5)÷6, (m+5)÷6, (n+5)÷6]
            end 


            ######################### three ############
            three=0 
            e=(l-1)÷3 
            f=(m-1)÷3 
            g=(n-1)÷3 
            
            S= floor(Int, e) 
            Sp= floor(Int, f) 
            Spp= floor(Int, g) 

            for s=0:S 
                for sp=0:Sp 
                    for spp=0:Spp  

                        if (l+2-(3*s))%3==0 && (m+2-(3*sp))%3==0 && (n+2-(3*spp))%3==0  
                            three += a[s+1,sp+1,spp+1] * a[(l+2-(3*s))÷3, (m-(3*sp)+2)÷3, (n-(3*spp)+2)÷3 ] 
                        end 
                        

                    end 
                end 
            end 

            
            ######################### four ############
            four=0 
            e=(l-1)÷2 
            f=(m-1)÷2 
            g=(n-1)÷2 
            
            S= floor(Int, e) 
            Sp= floor(Int, f) 
            Spp= floor(Int, g)  

            for s=0:S 
                for sp=0:Sp 
                    for spp=0:Spp 

                        ep=(l-1-(2*s))÷2 
                        fp=(m-1-(2*sp))÷2 
                        gp=(n-1-(2*spp))÷2 
                        
                        T= floor(Int, ep) 
                        Tp= floor(Int, fp) 
                        Tpp= floor(Int, gp) 

                        for t=0:T 
                            for tp=0:Tp 
                                for tpp=0:Tpp 

                                    if (l+1-(2*s)-(2*t))%2==0 && (m+1-(2*sp)-(2*tp))%2==0 && (n+1-(2*spp)-(2*tpp))%2==0 

                                        four += a[s+1,sp+1,spp+1] * a[t+1,tp+1,tpp+1] * 
                                        a[(l+1-(2*s)-(2*t))÷2, (m+1-(2*sp)-(2*tp))÷2, (n+1-(2*spp)-(2*tpp))÷2 ]  

                                    end 


                                end 
                            end 
                        end 



                    end 
                end 
            end 


            ######################### five ############ 
            five=0 
            
            S= l-1
            Sp= m-1
            Spp= n-1 

            for s=0:S 
                for sp=0:Sp 
                    for spp=0:Spp  

                        ep=(l-1-s)÷2 
                        fp=(m-1-sp)÷2 
                        gp=(n-1-spp)÷2 
                        
                        T= floor(Int, ep)
                        Tp= floor(Int, fp) 
                        Tpp= floor(Int, gp)   

                        for t=0:T 
                            for tp=0:Tp 
                                for tpp=0:Tpp 

                                    epp= (l-1-s-(2*t))÷2 
                                    fpp= (m-1-sp-(2*tp))÷2 
                                    gpp= (n-1-spp-(2*tpp))÷2 
                                    
                                    U= floor(Int, epp) 
                                    Up= floor(Int, fpp) 
                                    Upp= floor(Int, gpp) 

                                    for u=0:U 
                                        for up=0:Up 
                                            for upp=0:Upp  

                                                five+= a[s+1,sp+1,spp+1] * a[t+1,tp+1,tpp+1] * a[u+1,up+1,upp+1] *
                                                a[(l-s-(2*t)-(2*u)), (m-sp-(2*tp)-(2*up)), (n-spp-(2*tpp)-(2*upp)) ]  

                                            end 
                                        end 
                                    end 


                                end 
                            end 
                        end 


                    end 
                end 
            end 



            ########################### Benzen Coeeficients ################
            d[l,m,n]= ( one + (2*two) + (2*three) + (4*four) + (3*five)) ÷ 12 
            if d[l,m,n]!=0 
                println( l+5,",",m-1,",", n-1, ":  ",  d[l,m,n] )
            end

            

        end 
    end 
end  








