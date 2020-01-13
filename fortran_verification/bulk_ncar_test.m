clear

altu=10;
altt=2;
altq=2;
imx=5;
for i=1%1:imx
us(i,1)=-24.0+(i+1)*20.0
   vs(i,1)=-26.0+(i+1)*20.0
   sat(i,1)=-20.0+(i+1)*20
   qar(i,1)=0.0+(i+1)*0.1
   slp(i,1)=(90183.0+800.0*(i+1))/100.0
   sst(i,1)=-20.0+(i+1)*16.0
   wdv(i,1)=(us(i,1)^2+vs(i,1)^2)^0.5
end


[wsx,wsy,qla,qsn,evp,tu,qu,dtu,dqu,w10n]=bulk_ncar(us,vs,sat,qar,slp,sst,altu,altt,altq)

