      parameter (nfd=600,ned=600,nmd=6,nbd=50,nxy=350,ndi=72,

     &            nncn=16,nnbr=30)

c      parameter (nfd=200,ned=200,nmd=6,nbd=50,nxy=200,ndi=72,

c    &            nncn=20,nnbr=40)

c

       common  aaa,bbb,nonc,nhar,mmu,nhh,nsum,itr,ntotal(2,nmd),

     &         npo1,npo2,npo3,npor,nom1,nom2,nom,ncn,nbr,

     &         nelem,nb,nbb,nre,ncc,

     &         df,dmu0,domeg,dpi,

     &         dh(nfd,nbd),dhc(nfd,nncn+nnbr),dhl(nncn+nnbr,nfd),

     &         dhcl(nncn+nnbr,nncn+nnbr),dk(nfd+nncn+nnbr),

     &         dkk(nmd,nfd+nncn+nnbr),da(nfd+nncn+nnbr),

     &         daa(nmd,nfd+nncn+nnbr),daaa(nmd,nfd+nncn+nnbr),

     &         dnds(nnbr,3),ntabl(nnbr,5),dtabl(nnbr,5),

     &         dc(ned),dg(ned),dv(nnbr,nmd,2),

     &         dzz(nnbr,nnbr),dyy(nnbr,nnbr),das(nncn,nnbr),

     &         ds(ned,7),dcs(ned),dbx(ned,nmd,2),dby(ned,nmd,2),

     &         nod(ned,4),xy(nxy,2),dcur(100,nmd,2),dre(100),

     &         nodr(100,2),dd(ned,nmd,2,2),dn(nmd,2,2),

     &         dreluc(0:22,0:1),dsine(0:22,0:ndi),dcosi(0:22,0:ndi)