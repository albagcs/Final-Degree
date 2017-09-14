function [outlierst1_ut,outlierst1_br1,outlierst1_br2,outlierst1_lv1,outlierst1_lv2,outlierst1_pl]=outliers_s1(detJgw,plmask83,lv1mask83,lv2mask83,br1mask83,br2mask83,uterusB3,tlpl,thpl,tllv,thlv,tlbr,thbr) 
%for twin data:  detJgw is the determinant of jacobian matrix;
%lv1mask,lv2mask,br1mask,br2mask,uterusmask are the masks for livers,
%brains and uterus;tlpl,thpl,tllv,thlv,tlbr,thbr are the thresholds for
%over compression and over expansion for placenta, liver, and brain 
%outlierst1_* are the volumes after outlier rejection
%%%function [gwuter,gwuter_0,gwbr1,gwbr1_0,gwbr2,gwbr2_0,gwlv1,gwlv1_0,gwlv2,gwlv2_0,gwpl,gwpl_0]=outliers_step1(detJgw,plmask,lv1mask,lv2mask,br1mask,br2mask,uterusmask,tlpl,thpl,tllv,thlv,tlbr,thbr)
[x,y,z,t]=size(detJgw);
edjgw_0=zeros(t,z);for a=1:t for d=1:z count0=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(uterusB3(b,c,d))<0  count0=count0+1;end,end,edjgw_0(a,d)=count0;end,end,end
edjgw_tlpl=zeros(t,z);for a=1:t for d=1:z counttlpl=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(uterusB3(b,c,d))<tlpl && double(detJgw(b,c,d,a)).*double(uterusB3(b,c,d))>0 counttlpl=counttlpl+1;end,end,edjgw_tlpl(a,d)=counttlpl;end,end,end
edjgw_thpl=zeros(t,z);for a=1:t for d=1:z countthpl=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(uterusB3(b,c,d))>thpl countthpl=countthpl+1;end,end,edjgw_thpl(a,d)=countthpl;end,end,end
for a=1:t gwuter(a)=sum(edjgw_thpl(a,:)+edjgw_tlpl(a,:)+edjgw_0(a,:));end
%for a=1:t gwuter_0(a)=sum(edjgw_0(a,:));end
outlierst1_ut=find(gwuter==0);

br1s(:,:,:)=smooth3(br1mask83(:,:,:),'box',7); %evaluate a bigger region
for b=1:x for c=1:y for d=1:z if br1s(b,c,d)>0 br1cs(b,c,d)=1; else br1cs(b,c,d)=0; end,end,end,end
edjgwbr1_0=zeros(t,z);for a=1:t for d=1:z count0=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(br1cs(b,c,d))<0  count0=count0+1;end,end,edjgwbr1_0(a,d)=count0;end,end,end
edjgwbr1_tlbr=zeros(t,z);for a=1:t for d=1:z counttlbr=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(br1cs(b,c,d))<tlbr && double(detJgw(b,c,d,a)).*double(br1cs(b,c,d))>0 counttlbr=counttlbr+1;end,end,edjgwbr1_tlbr(a,d)=counttlbr;end,end,end
edjgwbr1_thbr=zeros(t,z);for a=1:t for d=1:z countthbr=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(br1cs(b,c,d))>thbr countthbr=countthbr+1;end,end,edjgwbr1_thbr(a,d)=countthbr;end,end,end
for a=1:t gwbr1(a)=sum(edjgwbr1_thbr(a,:)+edjgwbr1_tlbr(a,:)+edjgwbr1_0(a,:));end
%for a=1:t gwbr1_0(a)=sum(edjgwbr1_0(a,:));end
outlierst1_br1=find(gwbr1==0);

br2s(:,:,:)=smooth3(br2mask83(:,:,:),'box',7);
for b=1:x for c=1:y for d=1:z if br2s(b,c,d)>0 br2cs(b,c,d)=1; else br2cs(b,c,d)=0; end,end,end,end
edjgwbr2_0=zeros(t,z);for a=1:t for d=1:z count0=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(br2cs(b,c,d))<0  count0=count0+1;end,end,edjgwbr2_0(a,d)=count0;end,end,end
edjgwbr2_tlbr=zeros(t,z);for a=1:t for d=1:z counttlbr=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(br2cs(b,c,d))<tlbr && double(detJgw(b,c,d,a)).*double(br2cs(b,c,d))>0 counttlbr=counttlbr+1;end,end,edjgwbr2_tlbr(a,d)=counttlbr;end,end,end
edjgwbr2_thbr=zeros(t,z);for a=1:t for d=1:z countthbr=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(br2cs(b,c,d))>thbr countthbr=countthbr+1;end,end,edjgwbr2_thbr(a,d)=countthbr;end,end,end
for a=1:t gwbr2(a)=sum(edjgwbr2_thbr(a,:)+edjgwbr2_tlbr(a,:)+edjgwbr2_0(a,:));end
%for a=1:t gwbr2_0(a)=sum(edjgwbr2_0(a,:));end
outlierst1_br2=find(gwbr2==0);

lv1s(:,:,:)=smooth3(lv1mask83(:,:,:),'box',9);
for b=1:x for c=1:y for d=1:z if lv1s(b,c,d)>0 lv1cs(b,c,d)=1; else lv1cs(b,c,d)=0; end,end,end,end
edjgwlv1_0=zeros(t,z);for a=1:t for d=1:z count0=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(lv1cs(b,c,d))<0  count0=count0+1;end,end,edjgwlv1_0(a,d)=count0;end,end,end
edjgwlv1_tllv=zeros(t,z);for a=1:t for d=1:z counttllv=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(lv1cs(b,c,d))<tllv && double(detJgw(b,c,d,a)).*double(lv1cs(b,c,d))>0 counttllv=counttllv+1;end,end,edjgwlv1_tllv(a,d)=counttllv;end,end,end
edjgwlv1_thlv=zeros(t,z);for a=1:t for d=1:z countthlv=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(lv1cs(b,c,d))>thlv countthlv=countthlv+1;end,end,edjgwlv1_thlv(a,d)=countthlv;end,end,end
for a=1:t gwlv1(a)=sum(edjgwlv1_thlv(a,:)+edjgwlv1_tllv(a,:)+edjgwlv1_0(a,:));end
%for a=1:t gwlv1_0(a)=sum(edjgwlv1_0(a,:));end
outlierst1_lv1=find(gwlv1==0);

lv2s(:,:,:)=smooth3(lv2mask83(:,:,:),'box',9);
for b=1:x for c=1:y for d=1:z if lv2s(b,c,d)>0 lv2cs(b,c,d)=1; else lv2cs(b,c,d)=0; end,end,end,end
edjgwlv2_0=zeros(t,z);for a=1:t for d=1:z count0=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(lv2cs(b,c,d))<0  count0=count0+1;end,end,edjgwlv2_0(a,d)=count0;end,end,end
edjgwlv2_tllv=zeros(t,z);for a=1:t for d=1:z counttllv=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(lv2cs(b,c,d))<tllv && double(detJgw(b,c,d,a)).*double(lv2cs(b,c,d))>0 counttllv=counttllv+1;end,end,edjgwlv2_tllv(a,d)=counttllv;end,end,end
edjgwlv2_thlv=zeros(t,z);for a=1:t for d=1:z countthlv=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(lv2cs(b,c,d))>thlv countthlv=countthlv+1;end,end,edjgwlv2_thlv(a,d)=countthlv;end,end,end
for a=1:t gwlv2(a)=sum(edjgwlv2_thlv(a,:)+edjgwlv2_tllv(a,:)+edjgwlv2_0(a,:));end
%for a=1:t gwlv2_0(a)=sum(edjgwlv2_0(a,:));end
outlierst1_lv2=find(gwlv2==0);

pls(:,:,:)=smooth3(plmask83(:,:,:),'box',9);
for b=1:x for c=1:y for d=1:z if pls(b,c,d)>0 plcs(b,c,d)=1; else plcs(b,c,d)=0; end,end,end,end
edjgwpl_0=zeros(t,z);for a=1:t for d=1:z count0=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(plcs(b,c,d))<0  count0=count0+1;end,end,edjgwpl_0(a,d)=count0;end,end,end
edjgwpl_tlpl=zeros(t,z);for a=1:t for d=1:z counttlpl=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(plcs(b,c,d))<tlpl && double(detJgw(b,c,d,a)).*double(plcs(b,c,d))>0 counttlpl=counttlpl+1;end,end,edjgwpl_tlpl(a,d)=counttlpl;end,end,end
edjgwpl_thpl=zeros(t,z);for a=1:t for d=1:z countthpl=0; for b=1:x for c=1:y  if double(detJgw(b,c,d,a)).*double(plcs(b,c,d))>thpl countthpl=countthpl+1;end,end,edjgwpl_thpl(a,d)=countthpl;end,end,end
for a=1:t gwpl(a)=sum(edjgwpl_thpl(a,:)+edjgwpl_tlpl(a,:)+edjgwpl_0(a,:));end
%for a=1:t gwpl_0(a)=sum(edjgwpl_0(a,:));end
outlierst1_pl=find(gwpl==0);
