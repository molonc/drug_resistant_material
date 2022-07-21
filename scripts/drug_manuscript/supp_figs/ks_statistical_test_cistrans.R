library(ggplot2)

p <- ggplot(mtcars, aes(mpg, wt)) +
  geom_point() +
  facet_grid(. ~ cyl) +
  theme(panel.spacing = unit(1, "lines"))
p


dat_text <- data.frame(
  label = c("4 cylinders", "6 cylinders", "8 cylinders"),
  cyl   = c(4, 6, 8)
)
p1 <- p + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = 1,
  vjust   = -2
)
p1
hjust <- 0.5
vjust <- 0.5
td <- expand.grid(
  hjust=c(0, 0.5, 1),
  vjust=c(0, 0.5, 1)
)
hj <- 0.5
vj <- 1
for(hj in td$hjust){
  for(vj in td$vjust){
    print(hj)
    print(vj)
    p1 <- p + geom_text(data=dat_text, aes(x = -Inf, y = 4, 
                                           label=label,hjust=hj, vjust=vj)) 
    p1
    
  }
}
p1 <- p + geom_text(data=dat_text, aes(x = -Inf, y = -Inf, label=label, angle=20, hjust=hjust, vjust=vjust)) 
p1

ToothGrowth$dose <- as.factor(ToothGrowth$dose)
# Change box plot colors by groups
ggplot(ToothGrowth, aes(x=dose, y=len, fill=supp)) +
  geom_boxplot()
# Change the position
p<-ggplot(ToothGrowth, aes(x=dose, y=len, fill=supp)) +
  geom_boxplot(position=position_dodge(1))
p
unique(ToothGrowth$supp)
dat_text <- data.frame(
  label = c("4 cylinders", "6 cylinders", "8 cylinders"),
  ds   = c(1, 2, 3),
  supp=c('VC','VC','VC')  
)
View(ToothGrowth)
  
p1 <- p + annotate("text", x = c(0.5, 1, 2), y = 35,
          label = c("label 1", "label 2","label 3") , color="orange",
          size=7 , angle=45, fontface="bold")
p1
p1 <- p + geom_text(data=dat_text, aes(x = ds, y = 35, label=label,
                                       hjust=0.5)) #, hjust=hjust, vjust=vjust
p1
dat_text <- data.frame(
  label = c("4 cylinders", "6 cylinders", "8 cylinders"),
  cyl   = c(4, 6, 8),
  x     = c(20, 27.5, 25),
  y     = c(4, 4, 4.5)
)

p1 <- p + geom_text(
  data    = dat_text,
  mapping = aes(x = x, y = y, label = label)
)
p1



library(dplyr)
library(ggplot2)
ToothGrowth <- ToothGrowth %>%
  dplyr::mutate(dose_desc=case_when(
    dose == 0.5 ~ "aaa",
    dose == 1 ~ "bbb",
    dose == 2 ~ "ccc",
    TRUE ~ 'other'
  ))
ToothGrowth$dose_desc <- as.factor(ToothGrowth$dose_desc)
t1 <- ToothGrowth
t1$series <- 'Pt1'
t2 <- ToothGrowth
t2$series <- 'Pt2'
t3 <- ToothGrowth
t3$series <- 'Pt3'


t <- dplyr::bind_rows(as.data.frame(t1),as.data.frame(t2),as.data.frame(t3))

# Change box plot colors by groups
# Change the position
class(t)
dim(t)

p<-ggplot(t, aes(x=dose_desc, y=len, fill=supp)) +
  geom_boxplot(position=position_dodge(1)) + 
  facet_grid(. ~ series) +
  theme(panel.spacing = unit(1, "lines"))
p

# p<-ggplot(ToothGrowth, aes(x=dose, y=len, fill=supp)) +
#   geom_boxplot(position=position_dodge(1))
# p

dat_text1 <- data.frame(
  label = c("*", " ","**"),
  ds   = c('aaa', 'bbb', 'ccc'),
  supp=c('VC','VC','VC')  
)

dat_text <- data.frame(
  label = c("*", " ","**"),
  ds   = c('aaa', 'bbb', 'ccc'),
  supp=c('VC','VC','VC')
  # series=c('Pt1','Pt2','Pt3')
)

dat_text2 <- data.frame(
  label = c("***", "*","*"),
  ds   = c('aaa', 'bbb', 'ccc'),
  supp=c('VC','VC','VC')
  # series=c('Pt1','Pt2','Pt3')
)

dat_text3 <- data.frame(
  label = c("**", "***",""),
  ds   = c('aaa', 'bbb', 'ccc'),
  supp=c('VC','VC','VC')
  # series=c('Pt1','Pt2','Pt3')
)

# dat_text1 <- dat_text
dat_text1$series <- 'Pt1'
# dat_text2 <- dat_text2
dat_text2$series <- 'Pt2'
# dat_text3 <- dat_text3
dat_text3$series <- 'Pt3'
anno <- dplyr::bind_rows(dat_text1,dat_text2,dat_text3)
# View(ToothGrowth)

# p1 <- p + annotate("text", x = c(1, 2, 3), y = 35,
#                    label = c("*", " ","**") , color="black",
#                    size=7 , fontface="bold") #angle=45, 
# p1
# p1 <- p + geom_text(data=dat_text, aes(x = ds, y = 35, label=label,
#                                        hjust=0.5)) #, hjust=hjust, vjust=vjust
# p1

p1 <- p + geom_text(data=anno, aes(x = ds, y = 35, label=label,
                                       hjust=0.5)) #, hjust=hjust, vjust=vjust
p1


p.values <- c(9.5e-15, 0.02)
signif <- stats::symnum(p.values, corr = FALSE, na = FALSE, 
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
