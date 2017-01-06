clear all;
close all;

LAI=load('TLAI_out.txt');
GPP=load('GPP_out.txt');
parms=load('parms_out.txt');
obs=load('US-NR1_GPP.txt');


good=find(LAI(6,1:1000) > 0);
alive=find(LAI(6,1:1000) > 3.0 & LAI(6,1:1000) < 5.0);

plot(LAI(6,:),GPP(6,:),'rx');
for y=1:10
  for m=1:12
    GPP_av(m)=mean(mean(GPP(m:12:120,alive)));
    GPP_av_good(m)=mean(mean(GPP(m:12:120,good)));
    GPP_av_def(m)=mean(mean(GPP(m:12:120,1)));
    for n=1:1000
      GPP_seas_av(n,m)=mean(GPP(m:12:120,n));
    end;
    GPP_seas_ava(1:length(alive),m)=sort(GPP_seas_av(alive,m));
    GPP_u95a(m)=GPP_seas_ava(fix(length(alive)*0.975),m);
    GPP_l95a(m)=GPP_seas_ava(fix(length(alive)*0.025),m);

    GPP_seas_avg(1:length(good),m)=sort(GPP_seas_av(good,m));
    GPP_u95g(m)=GPP_seas_avg(fix(length(good)*0.975),m);
    GPP_l95g(m)=GPP_seas_avg(fix(length(good)*0.025),m);
  
    GPP_obs(m)=mean(obs(12+m:12:120,3));
    GPP_obs_unc(m)=mean(obs(12+m:12:120,5));
    for n=1:23
      mindex=(12+m):12:120;
      mygood=find(obs(mindex,5+n) > -900);
      if (length(mygood) > 0)
        GPP_obs_mod(m,n)=mean(obs(mindex(mygood),5+n));
      end
end;
end;
end;

plot(0.5:11.5,GPP_av*86400,'LineWidth',3)
hold on;
jbfill(0.5:11.5,(GPP_u95a)*86400,(GPP_l95a)*86400,'y','k',0.03)
hold on;
%jbfill(0.5:11.5,(GPP_u95g)*86400,(GPP_l95g)*86400,'c','k',0.1)
%       hold on;
%plot(GPP_av_good*86400,'r')
  plot(0.5:11.5, GPP_av_def*86400,'r','LineWidth',3)
  errorbar(0.5:11.5,GPP_obs/30.5,(-1*GPP_obs_unc*1.96)/30.5,(GPP_obs_unc*1.96)/30.5,'k','LineWidth',2.5)
  grey=[0.4,0.4,0.4];
plot(0.5:11.5,GPP_obs_mod(:,1:9)/30.5,'Color',grey)
  plot(0.5:11.5,GPP_obs_mod(:,11:22)/30.5,'Color',grey)
  legend('Ensemble mean','95% CI','Default CLM-CN','observations','NACP models')
  ah=gca;
  set(ah,'fontsize',14)
    axis([0 12 0 16])
  xlabel('Month')
  ylabel('GPP (gC m^-^2 d^-^1)')
