%example 
detJgw=zeros(110,110,80,82);
tlpl = 0.5; thpl=1.5; tllv=0.7; thlv=1.3; tlbr=0.8; thbr=1.2;
load_untouch_nii(sprintf('plmask83.nii'));plmask83=ans.img;
load_untouch_nii(sprintf('lv1mask83.nii'));lv1mask83=ans.img;
load_untouch_nii(sprintf('lv2mask83.nii'));lv2mask83=ans.img;
load_untouch_nii(sprintf('br1mask83.nii'));br1mask83=ans.img;
load_untouch_nii(sprintf('br2mask83.nii'));br2mask83=ans.img;
load_untouch_nii(sprintf('012115/uterus_B3.nii'));uterusB3=ans.img;

for a=1:82
    disp('Hola, mundo')
    load_nii(sprintf('elastix/oddevenB01/tmeosh_%d/spatialJacobian.nii',a));detJgw_B14(:,:,2:2:80,a)=ans.img;
    load_nii(sprintf('elastix/oddevenB01/tmoosh_%d/spatialJacobian.nii',a));detJgw_B14(:,:,1:2:79,a)=ans.img;
end
[outlierst1_ut,outlierst1_br1,outlierst1_br2,outlierst1_lv1,outlierst1_lv2,outlierst1_pl]=outliers_s1(detJgw_B14,plmask83,lv1mask83,lv2mask83,br1mask83,br2mask83,uterusB3,tlpl,thpl,tllv,thlv,tlbr,thbr);


