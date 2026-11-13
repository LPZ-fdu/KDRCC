# Run RUN_EXAMPLE.R first. This script produces ggplot2 figures from data frames.
script<-NULL
for(fr in sys.frames())if(!is.null(fr$ofile))script<-fr$ofile
args<-grep("^--file=",commandArgs(FALSE),value=TRUE)
if(is.null(script)&&length(args))script<-sub("^--file=","",args[1])
root<-if(is.null(script))getwd() else dirname(normalizePath(script,winslash="/",mustWork=TRUE))
library(ggplot2)
out<-file.path(root,"output")
single<-read.csv(file.path(out,"single_effects.csv"))
joint<-read.csv(file.path(out,"joint_effect.csv"))
contribution<-read.csv(file.path(out,"contribution_shares.csv"))
p1<-ggplot(single,aes(x=concentration,y=or))+
  geom_ribbon(aes(ymin=or_lower,ymax=or_upper),fill="#0072B2",alpha=.2)+
  geom_line(colour="#0072B2",size=.7)+geom_hline(yintercept=1,linetype=3)+
  facet_wrap(~exposure,scales="free")+labs(x="Exposure concentration",y="Odds ratio")+
  theme_classic(base_size=14)
p2<-ggplot(joint,aes(position,or))+
  geom_ribbon(aes(ymin=or_lower,ymax=or_upper),fill="#0072B2",alpha=.2)+
  geom_line(colour="#0072B2",size=.7)+geom_hline(yintercept=1,linetype=3)+
  labs(x="Common exposure quantile",y="Odds ratio")+theme_classic(base_size=14)
p3<-ggplot(contribution,aes(x=reorder(variable,share_percent),y=share_percent))+
  geom_col(fill="#0072B2")+coord_flip()+
  labs(x="Exposure",y="Mean absolute contribution (%)")+theme_classic(base_size=14)
ggsave(file.path(out,"custom_single.png"),p1,width=10,height=7,dpi=300)
ggsave(file.path(out,"custom_joint.png"),p2,width=8,height=6,dpi=300)
ggsave(file.path(out,"custom_contribution.png"),p3,width=8,height=6,dpi=300)
