############################################################################################
########################## Descriptive Analysis Recurrent Events ###########################
############################################################################################

## [1] Plot 1: recurrent events encoding example
## [2] Plot 2: multiple events & composite endpoint
## [3] Plot 3: Passing to Gap Times may cause loss of information
## [4] Plot 4: Gap times vs Number event
## [5] Plot 5: Gap Times Histogram

## Plot 1: recurrent events encoding example
require(gridExtra)
df_plot1<-data.frame(
  COD_REG=c("subject 1","subject 1","subject 1","subject 2","subject 3",
            "subject 3","subject 3","subject 3","subject 4"),
  event =c("event","event","censoring","censoring","event","event",
           "event","censoring","censoring"),
  time_y=c(0.2,0.3,0.7,2,0.35,0.7,1.3,2,1),
  left=rep(0.0,9),
  right=c(0.7,0.7,0.7,2.0,2.0,2.0,2.0,2.0,1.0)
)

df_plot2<-data.frame(
  COD_REG=c("subject 1","subject 2","subject 3","subject 4"),
  event =c("event","censoring","event","censoring"),
  time_y=c(0.2,2,0.35,1),
  left=c(0,0,0,0)
)

plot1 <-ggplot(data=df_plot1,aes(x=COD_REG,y=as.double(time_y),group=COD_REG)) +
  geom_segment(aes(x=COD_REG,y=left,xend=COD_REG,yend=right))+
  geom_point(aes(color=event),size=3)+
  scale_color_manual(values=c("white", "black"))+
  xlab("")+
  ylab("years")+
  coord_flip()+
  geom_hline(yintercept = 2,linetype= "dotted")+
  theme(legend.title = element_blank())

plot2 <- 
  ggplot(data=df_plot2,aes(x=COD_REG,y=as.double(time_y),group=COD_REG)) +
  geom_segment(aes(x=COD_REG,y=left,xend=COD_REG,yend=time_y))+
  geom_point(aes(color=event),size=3)+
  scale_color_manual(values=c("white", "black"))+
  xlab("")+
  ylab("years")+
  coord_flip()+
  geom_hline(yintercept = 2,linetype= "dotted")+
  theme(legend.title = element_blank())

grid.arrange(plot1, plot2, ncol=1)

rm(df_plot1,df_plot2,plot1,plot2)


## Plot 2: multiple events & composite endpoint
data_head=data[1:18,]
data_head$time_y=data_head$timeEvent/365
data_head$event=factor(data_head$event, labels=c('terminal','recurrent'))
data_head$left=rep(0.0,18)
ggplot(data=data_head,aes(x=COD_REG,y=as.double(time_y),group=COD_REG)) +
  geom_segment(aes(x=COD_REG,y=left,xend=COD_REG,yend=time_y))+
  geom_point(aes(color=event,shape=status),size=3)+
  xlab("ID")+
  ylab("Time [years]")+
  coord_flip()
rm(data_head)

## Plot 3: Passing to Gap Times may cause loss of information
data_plot3=data[COD_REG %in% c(10007000,10000717),]
data_plot3$time_y=data_plot3$GapEvent
data_plot3$event=factor(data_plot3$event, labels=c('terminal','recurrent'))
data_plot3$left=rep(0.0,11)
data_plot3$index=c(1.5,2.0,2.5,3.0,3.5,4.0,4.5,6.5,7.0,7.5,8.0)
ggplot(data=data_plot3,aes(x=index,y=as.double(time_y))) +
  geom_segment(aes(x=index,y=left,xend=index,yend=time_y))+
  geom_point(aes(color=event,shape=status),size=3)+
  xlab("patient")+
  ylab("Time [days]")+
  scale_x_discrete(breaks=NULL)+
  annotate("text",x=3.0,y=-20,label="10000717", size=3)+
  annotate("text",x=7.25,y=-20,label="10007000",size=3)
rm(data_plot3)

## Plot 4: Gap times vs Number event
library(plyr)
print(table(data$hosp<=9))
hist(data[hosp>9]$hosp)
indiv.stats <- ddply(data[hosp<=9,],c("COD_REG","hosp"),summarize,time2event = GapEvent)
indiv.stats  <- indiv.stats[order(indiv.stats$hosp, -indiv.stats$time2event),]
indiv.stats$order <- c(1:length(indiv.stats$COD_REG))
indiv.stats <- na.omit(indiv.stats)
medianepi <- ddply(indiv.stats,"hosp",summarize, medep = median(order))
medianepi <- na.omit(medianepi)
ggplot(data=indiv.stats, aes(x=time2event, y=order)) +
  geom_rect(mapping = aes(xmin=0,xmax=indiv.stats$time2event,ymin=indiv.stats$order,
                          ymax=indiv.stats$order+1,fill=factor(hosp))) +
  scale_y_continuous(breaks=NULL) +
  xlab("Time [days]") + ylab("") + labs(fill="Event\nNumber") + theme_gray()
rm(indiv.stats,medianepi)

## Plot 5: Gap Times Histogram
p2<-ggplot(dataDeath, aes(x=GapEvent)) +
  geom_histogram(position="identity", binwidth=30)+
  xlab("Terminal Gap time [days]")+
  ylab("count")

p1<-ggplot(data, aes(x=GapEvent)) +
  geom_histogram(position="identity", binwidth=30)+
  xlab("Gap time [days]")+
  ylab("count")

ggarrange(p1,p2,labels = c("",""),nrow=1)