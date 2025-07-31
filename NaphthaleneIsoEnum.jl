

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
a = count_acyclic_isomers(8,4,18); #typeof(a): Array{BigInt, 3} 


####################### Nephthalene ############################
####################### Nephthalene ############################ 
####################### Nephthalene ############################ 

#Enter L carbon, M oxygen, and N Hydrogen atom
L=3     # L_min = 1
M=3     # M_min = 1
N= 12   # N_min  =9

########## Term1_Naphthalene_coeffs #################################### 
first=zeros(BigInt,L,M,N); 
for l=1:L 
    for m=1:M 
        for n=1:N 

            one=0 

            S= l-1  
            Sp= m-1   
            Spp= n-1 

            for s=0:S 
                for sp=0:Sp 
                    for spp=0:Spp  

                        T=l-s-1 
                        Tp=m-sp-1 
                        Tpp=n-spp-1 

                        for t=0:T 
                            for tp=0:Tp 
                                for tpp=0:Tpp 

                                    U= l-s-t-1 
                                    Up= m-sp-tp-1 
                                    Upp= n-spp-tpp-1 

                                    for u=0:U 
                                        for up=0:Up 
                                            for upp=0:Upp 

                                                V= l-s-t-u-1
                                                Vp= m-sp-tp-up-1 
                                                Vpp= n-spp-tpp-upp-1 

                                                for v=0:V 
                                                    for vp=0:Vp 
                                                        for vpp=0:Vpp 

                                                            W= l-s-t-u-v-1 
                                                            Wp= m-sp-tp-up-vp-1 
                                                            Wpp= n-spp-tpp-upp-vpp-1 

                                                            for w=0:W 
                                                                for wp=0:Wp 
                                                                    for wpp=0:Wpp  

                                                                        X= l-s-t-u-v-w-1  
                                                                        Xp= m-sp-tp-up-vp-wp-1 
                                                                        Xpp= n-spp-tpp-upp-vpp-wpp-1 

                                                                        for x=0:X 
                                                                            for xp=0:Xp 
                                                                                for xpp=0:Xpp 

                                                                                    Y= l-s-t-u-v-w-x-1 
                                                                                    Yp= m-sp-tp-up-vp-wp-xp-1 
                                                                                    Ypp= n-spp-tpp-upp-vpp-wpp-xpp-1 

                                                                                    for y=0:Y 
                                                                                        for yp=0:Yp 
                                                                                            for ypp=0:Ypp  

                                                                                                one +=  a[s+1,sp+1,spp+1] * a[t+1,tp+1,tpp+1] * a[u+1,up+1,upp+1] * a[v+1,vp+1,vpp+1] * 
                                                                                                a[w+1,wp+1,wpp+1] * a[x+1,xp+1,xpp+1] *a[y+1,yp+1,ypp+1] * 
                                                                                                a[l-s-t-u-v-w-x-y,m-sp-tp-up-vp-wp-xp-yp,n-spp-tpp-upp-vpp-wpp-xpp-ypp] 

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


                                end 
                            end 
                        end 



                    end 
                end 
            end 
            
            
            first[l,m,n]= one 
            #if first[l,m,n]!=0 
                #println( l-1,",",m-1,",", n-1, ":    ", first[l,m,n] )
            #end


        end 
    end 
end  

########## Term2_Naphthalene_coeffs ####################################
second=zeros(BigInt,L,M,N); 

for l=1:L 
    for m=1:M 
        for n=1:N 

            two=0  

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
                        Tpp=floor(Int, gp) 

                        for t=0:T 
                            for tp=0:Tp 
                                for tpp=0:Tpp  

                                    epp= (l-1-(2*s)-(2*t))÷2    
                                    fpp= (m-1-(2*sp)-(2*tp))÷2   
                                    gpp= (n-1-(2*spp)-(2*tpp))÷2   

                                    U=floor(Int, epp)  
                                    Up=floor(Int, fpp)  
                                    Upp=floor(Int, gpp) 


                                    for u=0:U 
                                        for up=0:Up 
                                            for upp=0:Upp  

                                                if (l+1-(2*s)-(2*t)-(2*u))%2==0 && (m+1-(2*sp)-(2*tp)-(2*up))%2==0 && (n+1-(2*spp)-(2*tpp)-(2*upp))%2==0 
                                                    two+= a[s+1,sp+1,spp+1] * a[t+1,tp+1,tpp+1] * a[u+1,up+1,upp+1] * 
                                                    a[(l+1-(2*s)-(2*t)-(2*u))÷2 , (m+1-(2*sp)-(2*tp)-(2*up))÷2 , (n+1-(2*spp)-(2*tpp)-(2*upp))÷2 ] 

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


            second[l,m,n]= two 
            #if second[l,m,n]!=0 
                #println( l-1,",",m-1,",", n-1, ":    ", second[l,m,n] )
            #end
            

        end 
    end 
end 



############# Naphthalene_coeffs #####################

c=zeros(BigInt,L,M,N); 
for l=1:L 
    for m=1:M 
        for n=1:N 

            c[l,m,n]= (first[l,m,n] .+ (3*second[l,m,n])) .÷ 4 
            if c[l,m,n]!=0 
                println( l+9,",",m-1,",", n-1, ":    ",  c[l,m,n] )
            end

             
        end 
    end 
end 




