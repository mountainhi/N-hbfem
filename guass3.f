c ***************************************************************       
c ********  gaussian elemination method, gauss  ver.3.2 *********       
c                                                                       
      subroutine gauss3                                                 
c                                                                       
       include 'common.inc'
c                                                                       
c  number of matrix width                                               
      ndia=nb                                                           

c                                                                       
c *****************************************************************     

c     dk(nom)        : forced vector                                    

c     dhc(nom1,nom2) : side matrix of dh matrix                         

c     dhl(nom2,nom1) : under matrix of dh matrix                        

c     dhcl(nom2,nom2): coner matrix of dh matrix                        

c     daa(nom)       : potential      (nob)                             

c     nom            : number of unknown potential                      

c              nom=nom1+nom2                                            

c     nom1           :                                                  

c     nom2           :                                                  

c     ndia           :                                                  

c              ndia=nb                                                  

c     nbb            : band width                                       

c              nbb=2*nb-1                                               

c                                                                       

c ************** eliminate l of matrix dh(i,j) ********************     

c                                                                       

              nbb=2*nb-1                                                

c                                                                       

      do 100 nn=1, nom-1                                                

c                                                                       

        if (nn.gt.nom1)                      go to 110                  

        dt1=1.0/dh(nn,ndia)                                             

c ***** dh ******                                                       

c (nm,nn)                                                               

        do 200 nm=nn+1, nom                                             

          if(nm.gt.nom1)                     go to 210                  

c ***** dh-dh ******                                                    

          j1=ndia-nm+nn                                                 

          i1=nm                                                         

          if(j1.lt.1)                        go to 200                  

c                                                                       

          if(dh(i1,j1).eq.0.0)               go to 200                  

          dt=dt1*dh(i1,j1)                                              

c                                                                       

          dk(i1)=dk(i1)-dt*dk(nn)                                       

c (nm,nk)                                                               

          ii1=nm                                                        

          ii2=nn                                                        

          do 300 nk=nn, nom                                             

c                                                                       

            if(nk.gt.nom1)                   go to 310                  

c  ****** dh-dh *****                                                   

            jj1=ndia-nm+nk                                              

            jj2=ndia-nn+nk                                              

            if((jj1.gt.nbb).or.(jj1.lt.1))    go to 300                 

            if((jj2.gt.nbb).or.(jj2.lt.1))    go to 300                 

c                                                                       

            dh(ii1,jj1)=dh(ii1,jj1)-dt*dh(ii2,jj2)                      

            go to 300                                                   

c   ***** hdc-hdc *****                                                 

  310       continue                                                    

            jj1=nk-nom1                                                 

            jj2=jj1                                                     

            dhc(ii1,jj1)=dhc(ii1,jj1)-dt*dhc(ii2,jj2)                   

c      -- nk --                                                         

  300     continue                                                      

c                                                                       

        go to 200                                                       

c                                                                       

c   ***** dh-dhl, dhc-dhlc *****                                        

  210     continue                                                      

            ii1=nm-nom1                                                 

            ii2=nn                                                      

            jj1=nn                                                      

c                                                                       

            if(dhl(ii1,jj1).eq.0.0)           go to 200                 

            dt=dt1*dhl(ii1,jj1)                                         

c                                                                       

            i1=ii1+nom1                                                 

            dk(i1)=dk(i1)-dt*dk(ii2)                                    

c (nm,nk)                                                               

          do 400 nk=nn, nom                                             

            if(nk.gt.nom1)                    go to 410                 

c  ****** dh-dhl *****                                                  

            jj1=nk                                                      

            jj2=ndia-nn+nk                                              

            if(jj2.gt.nbb)                    go to 400                 

c                                                                       

            if(dh(ii2,jj2).eq.0.0)            go to 400                 

            dhl(ii1,jj1)=dhl(ii1,jj1)-dt*dh(ii2,jj2)                    

            go to 400                                                   

c   ***** hdc-hdcl *****                                                

  410     continue                                                      

            jj1=nk-nom1                                                 

            jj2=jj1                                                     

            if(dhc(ii2,jj2).eq.0.0)           go to 400                 

            dhcl(ii1,jj1)=dhcl(ii1,jj1)-dt*dhc(ii2,jj2)                 

c      -- nk --                                                         

  400     continue                                                      

c      -- nm --                                                         

  200   continue                                                        

c                                                                       

        go to 100                                                       

c                                                                       

c ****** dhlc *******                                                   

  110   continue                                                        

c (nm,nn)                                                               

        nnn=nn-nom1                                                     

c ***** pipotting *****                                                 

        l=nnn                                                           

        al=abs(dhcl(nnn,nnn))                                           

        do 20 j=nnn+1, nom2                                             

          if(abs(dhcl(j,nnn)).gt.al) then                               

            l=j                                                         

            al=abs(dhcl(j,nnn))                                         

          end if                                                        

   20   continue                                                        

c                                                                       

        if (nnn.ne.l) then                                              

          do 30 k=nnn,nom2                                              

            dtmp=dhcl(nnn,k)                                            

            dhcl(nnn,k)=dhcl(l,k)                                       

            dhcl(l,k)=dtmp                                              

   30     continue                                                      

          dtmp=dk(nom1+nnn)                                             

          dk(nom1+nnn)=dk(nom1+l)                                       

          dk(nom1+l)=dtmp                                               

        end if                                                          

c       print '( i4, e14.4)',nnn,dhcl(nnn,nnn)                          

c       if(abs(dhcl(nnn,nnn)).gt.1.0e-20)  go to 1110                   

c       print '(''bunnbo zero'',i4)',nnn                                

 1110   continue                                                        

c                                                                       

        dt1=1.0/dhcl(nnn,nnn)                                           

        do 500 nm=nn+1,nom                                              

c                                                                       

          jj2=nn-nom1                                                   

          ii2=nm-nom1                                                   

          if(dhcl(ii2,jj2).eq.0.0)                go to 500             

          dt=dt1*dhcl(ii2,jj2)                                          

c                                                                       

          ii1=nm                                                        

          dk(ii1)=dk(ii1)-dt*dk(nn)                                     

c (nm,nk)                                                               

            ii1=nm-nom1                                                 

            ii2=nn-nom1                                                 

          do 600 nk=nn, nom                                             

            jj1=nk-nom1                                                 

            jj2=jj1                                                     

c                                                                       

            dhcl(ii1,jj1)=dhcl(ii1,jj1)-dt*dhcl(ii2,jj2)                

c      -- nk --                                                         

  600     continue                                                      

c      -- nm --                                                         

  500   continue                                                        

c      -- nn --                                                         

  100 continue                                                          

c                                                                       

c ******************* calculate da(i) *******************               

c                                                                       

      da(nom)=dk(nom)/dhcl(nom2,nom2)                                   

c                                                                       

      do 800 nn=1, nom-1                                                

        ii1=nom-nn                                                      

        if(ii1.le.nom1)                          go to 810              

c                                                                       

        ii2=ii1-nom1                                                    

        do 700 nk=ii2+1, nom2                                           

          jj2=nk+nom1                                                   

          dk(ii1)=dk(ii1)-dhcl(ii2,nk)*da(jj2)                          

  700 continue                                                          

        da(ii1)=dk(ii1)/dhcl(ii2,ii2)                                   

        go to 800                                                       

c                                                                       

  810   continue                                                        

        do 900 nk=ii1+1, nom                                            

          if(nk.gt.nom1)                         go to 910              

c  ***** dh *****                                                       

          jj1=ndia-ii1+nk                                               

          if(jj1.gt.nbb)                          go to 900             

          dk(ii1)=dk(ii1)-dh(ii1,jj1)*da(nk)                            

          go to 900                                                     

c                                                                       

  910   continue                                                        

          jj1=nk-nom1                                                   

          dk(ii1)=dk(ii1)-dhc(ii1,jj1)*da(nk)                           

c   -- nk --                                                            

  900 continue                                                          

         if(dh(ii1,ndia).eq.0.0)       then                             

            print '('' matrix error (3)'',i5)',nn                       

         end if                                                         

         da(ii1)=dk(ii1)/dh(ii1,ndia)                                   

c   -- nn --                                                            

  800 continue                                                          

       return                                                           

       end     
                                                         
                                                 