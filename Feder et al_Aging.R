## LOAD REQUIRED PACKAGES
library(dbplyr)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(reshape2)
library(rstanarm)
library(viridis)
library(patchwork)
library(brms)
theme_set(theme_classic(base_size = 16))

## LOAD DATA
setwd("~/Desktop/SPRF/Geladas/Aging/GitHub")
life.table<-read.csv("Life_Tables.csv")
ibi.chacma<-read.csv("Chacma_IBI.csv")
ibi.gelada<-read.csv("Gelada_IBI.csv")
surv.chacma<-read.csv("Chacma_Surv.csv")
surv.gelada<-read.csv("Gelada_Surv.csv")

## PART 1: LIFE TABLE ANALYSES
## PLOTS
color.palette<-c("#440154FF", "#21908CFF", "darkgray", "grey40")

plot.fertility<-ggplot(data=life.table, aes(x=Age,y=mx, group=Population, color=Population)) + 
  scale_color_manual(values=color.palette) + ylab("Fertility\n(births/female/year)") + xlab("Age (years)") +
  theme(legend.position = "none") +
  geom_line(lwd=1, aes(linetype=Population)) + scale_linetype_manual(values=c("Gelada - Simiens"="solid", "Chacma - Moremi"="solid","Yellow - Amboseli"="dotted", "Olive - Gombe"="dotted"))

plot.survival<-ggplot(data=life.table, aes(x=Age,y=lx, group=Population, color=Population)) + 
  scale_color_manual(values=color.palette) + ylab("Survivorship") + xlab("Age (years)") + 
  geom_line(lwd=1, aes(linetype=Population)) + scale_linetype_manual(values=c("Gelada - Simiens"="solid", "Chacma - Moremi"="solid","Yellow - Amboseli"="dotted", "Olive - Gombe"="dotted")) +
  theme(legend.position =  c(0.65,0.85), legend.text = element_text(size=12), legend.key.width = unit(1,"cm"))

plot.fertility + plot.survival + plot_annotation(tag_levels = "a")
ggsave("Fertility_and_survival_chacma_v_geladas.jpg", units='cm', height=10, width=24)

## CALCULATE EXPECTED OFFSPRING
sum(life.table$mx[life.table$Population=="Chacma - Moremi"]*life.table$lx[life.table$Population=="Chacma - Moremi"]/0.65)
sum(life.table$mx[life.table$Population=="Gelada - Simiens"]*life.table$lx[life.table$Population=="Gelada - Simiens"])

## CALCULATE EXPECTED OFFSPRING FOR THOSE THAT SURVIVE INFANCY
sum(life.table$mx[life.table$Population=="Chacma - Moremi"]*life.table$lx[life.table$Population=="Chacma - Moremi"]/0.45)
sum(life.table$mx[life.table$Population=="Gelada - Simiens"]*life.table$lx[life.table$Population=="Gelada - Simiens"]/0.69)

## RUN MODELS
## FERTILITY
life.table$Births[is.na(life.table$Births)]<-round(life.table$mx[is.na(life.table$Births)]*life.table$Entered[is.na(life.table$Births)])

library(brms)
fertility.model<-brm(data=life.table, Births|trials(Entered)~poly(Age,2)*Population, 
            prior=c(
              prior(normal(0, 5), "Intercept"),
              prior(normal(0, 1), "b")),
            family="binomial", chains=4,iter=4000)
summary(fertility.model)
conditional_effects(fertility.model, "Age:Population")
conditional_effects(fertility.model, "Population")

## MORTALITY
mortality.model<-brm(data=life.table, Died|trials(Entered)~poly(Age,2)*Population, 
                     prior=c(
                       prior(normal(0, 5), "Intercept"),
                       prior(normal(0, 1), "b")),
                     family="binomial", chains=4,iter=4000)
summary(mortality.model)
conditional_effects(mortality.model, "Age:Population")
conditional_effects(mortality.model, "Population")

## PART 2: INTERBIRTH INTERVALS
pu.min<-min(ibi.chacma$ibi.num[ibi.chacma$ibi.censored==1])
tg.min<-min(ibi.gelada$ibi.num[ibi.gelada$ibi.censored==1])

## RUN MODELS
fit.pu <- stan_surv(Surv(ibi.adj, ibi.censored) ~ tve(mom.z) + tve(mom.z.2) + (1|mom.id) + (1|year), 
                    data=ibi.chacma, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.pu, prob=c(0.025,0.975)), digits=2)

fit.tg <- stan_surv(Surv(ibi.adj, ibi.censored) ~ tve(mom.z) + tve(mom.z.2) + (1|mom.id) + (1|year) + (1|unit), 
                    data=ibi.gelada, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.tg, prob=c(0.025,0.975)), digits=2)

## PARITY MODELS
fit.pu.parity <- stan_surv(Surv(ibi.adj, ibi.censored) ~ tve(mom.z) + parity + (1|mom.id) + (1|year), 
                    data=ibi.chacma, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.pu.parity, prob=c(0.025,0.975)), digits=2)

fit.tg.parity <- stan_surv(Surv(ibi.adj, ibi.censored) ~ tve(mom.z) + parity + (1|mom.id) + (1|year) + (1|unit), 
                    data=ibi.gelada, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.tg.parity, prob=c(0.025,0.975)), digits=2)

## COMPARE CHACMA MODELS
loo1<-loo(fit.pu)
loo2<-loo(fit.pu.parity)
loo_compare(loo1, loo2)

## COMPARE GELADA MODELS
loo1<-loo(fit.tg)
loo2<-loo(fit.tg.parity)
loo_compare(loo1, loo2)

## MAKE HR DATASET FOR CHACMA
try1<-plot(fit.pu, plotfun="tve",pars="mom.z")
try1<-try1$data
try2<-plot(fit.pu, plotfun="tve",pars="mom.z.2")
try2<-try2$data

colnames(try1)<-c("Time", "HR","LB", "UB")
colnames(try2)<-c("Time", "HR","LB", "UB")
try1$Covariate<-"Age.z"
try2$Covariate<-"Age.z.2"

hazards<-rbind(try1,try2)

## PLOT HAZARDS CHANGE
cols<-c("#1F78B4", "#33A02C")
hr1<-ggplot(data=hazards,aes(x=(Time+pu.min)/365.25, y=log(HR), group=Covariate, fill=Covariate)) +
  geom_vline(aes(xintercept=pu.min/365.25),lty=2) +
  annotate("rect", xmin=0, xmax=pu.min/365.25, ymin=-6, ymax=2, alpha=0.5, fill="lightgray")+
  scale_fill_manual(values=cols) + geom_hline(aes(yintercept=0), lty=2) +
  geom_ribbon(aes(ymin=log(LB), ymax=log(UB)), alpha=0.5) + geom_line(color="black") +
  scale_y_continuous(expand=c(0,0), breaks=seq(-6,2,2)) + scale_x_continuous(expand=c(0,0), limits=c(0,4),breaks = seq(0,4,2)) +
  xlab("years since previous birth") + ylab("log(Hazard ratio)") +
  ggtitle("Chacma baboons") + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold"), legend.position = "none")
hr1

## MAKE HR DATASET FOR GELADA
try1<-plot(fit.tg, plotfun="tve",pars="mom.z")
try1<-try1$data
try2<-plot(fit.tg, plotfun="tve",pars="mom.z.2")
try2<-try2$data

colnames(try1)<-c("Time", "HR","LB", "UB")
colnames(try2)<-c("Time", "HR","LB", "UB")
try1$Covariate<-"Age.z"
try2$Covariate<-"Age.z.2"

hazards<-rbind(try1,try2)

## PLOT HAZARDS CHANGE
RColorBrewer::brewer.pal(4, "Paired")
hr2<-ggplot(data=hazards,aes(x=(Time+tg.min)/365.25, y=log(HR), group=Covariate, fill=Covariate)) +
  geom_vline(aes(xintercept=tg.min/365.25),lty=2) +
  annotate("rect", xmin=0, xmax=tg.min/365.25, ymin=-6, ymax=2, alpha=0.5, fill="lightgray")+ geom_hline(aes(yintercept=0), lty=2) +
  geom_ribbon(aes(ymin=log(LB), ymax=log(UB)), alpha=0.5) + geom_line(color="black") +
  scale_y_continuous(expand=c(0,0), breaks=seq(-6,2,2)) + scale_x_continuous(expand=c(0,0), limits=c(0,6),breaks = seq(0,6,2)) +
  xlab("years since previous birth") + ylab("") + scale_fill_manual(labels=expression("Maternal age","Maternal age"^2), values=cols) +
  ggtitle("Geladas") + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold"), legend.text.align = 0)
hr2

hr1+hr2 + plot_annotation(tag_levels = "a")
ggsave("TVE_chacmas_v_geladas.jpg", units='cm', height=10, width=24)

## PLOT PROJECTIONS -- CHACMAS
vals<-round(quantile(ibi.chacma$mom.z, probs=c(.1,.3,.5,.7,.9)),2)

age.effect <- dplyr::select(ibi.chacma, mom.id, year) %>%
  tidyr::expand(mom.z=c(vals),
                nesting(mom.id, year)) %>%
  dplyr::mutate(id = 1:n(),
                age.cats = as.factor(mom.z),mom.z.2=c(mom.z^2)) %>%
  nest_legacy(-c(age.cats)) %>%
  dplyr::mutate(., preds = map(
    data,
    ~ posterior_survfit(
      fit.pu, prob=0.89,
      newdata = .x,
      times = 0,
      standardise = TRUE,
      extrapolate = TRUE,
      dynamic = TRUE
    )
  ))


age.effect1 <- unnest(age.effect, preds) %>%
  dplyr::select(-data)

?hcl.colors
blues<-RColorBrewer::brewer.pal(8,"Blues")
cb<-c("#FF7F00",  "#FDBF6F" ,"gray","#A6CEE3", "#1F78B4")
cb2<-c("#008837", "#A6DBA0", "gray", "#C2A5CF", "#7B3294")
cb3<-c("#008837", "#A6DBA0", "gray", "#F1B6DA", "#D01C8B")
blues<-blues[4:8]
plot.pu1<-ggplot(age.effect1,
                 aes(
                   x=(time+tg.min)/365.25 , y=1-median,
                   color=age.cats
                 )) + 
  geom_vline(aes(xintercept=pu.min/365.25),lty=2) +
  labs(fill="age category") + theme(legend.position = "none") +
  annotate("rect", xmin=0, xmax=pu.min/365.25, ymin=0, ymax=1, alpha=0.5, fill="lightgray")+
  geom_line(lwd=1.25,alpha=1)  + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold")) +
  scale_color_manual(values = cb) + xlab("years since previous birth") + ylab("") + ggtitle("Chacma baboons") + ylab("proportion closing IBI") +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.5)) + scale_x_continuous(expand=c(0,0), limits=c(0,4),breaks = seq(0,4,2))
plot.pu1

## PLOT PROJECTIONS -- GELADAS
vals<-round(quantile(ibi.gelada$mom.z, probs=c(.1,.3,.5,.7,.9)),2)

age.effect <- dplyr::select(ibi.gelada, mom.id, year,unit) %>%
  tidyr::expand(mom.z=c(vals),
                nesting(mom.id, year,unit)) %>%
  dplyr::mutate(id = 1:n(),
                age.cats = as.factor(mom.z),mom.z.2=c(mom.z^2)) %>%
  nest_legacy(-c(age.cats)) %>%
  dplyr::mutate(., preds = map(
    data,
    ~ posterior_survfit(
      fit.tg, prob=0.89,
      newdata = .x,
      times = 0,
      standardise = TRUE,
      extrapolate = TRUE,
      dynamic = TRUE
    )
  ))


age.effect2 <- unnest(age.effect, preds) %>%
  dplyr::select(-data)

plot.tg1<-ggplot(age.effect2,
                 aes(
                   x=(time+tg.min)/365.25 , y=1-median,
                   color=age.cats
                 )) + 
  geom_vline(aes(xintercept=tg.min/365.25),lty=2) +
  labs(fill="age category") + theme(legend.position = "none") +
  annotate("rect", xmin=0, xmax=tg.min/365.25, ymin=0, ymax=1, alpha=0.5, fill="lightgray")+
  geom_line(lwd=1.25,alpha=1)  + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold")) +
  scale_color_manual(values = cb) + xlab("years since previous birth") + ylab("") + ggtitle("Geladas") +
  scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.5)) + scale_x_continuous(expand=c(0,0), limits=c(0,6),breaks = seq(0,6,2))
plot.tg1

## PART 2B: IBI AND TAKEOVERS IN GELADAS
fit.tg.TO <- stan_surv(Surv(ibi.adj, ibi.censored) ~ tve(mom.z) + tve(mom.z.2) + takeover + (1|mom.id) + (1|unit) + (1|year), 
                    data=ibi.gelada, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.tg.TO, prob=c(0.025,0.975)), digits=2)

fit.tg.NO.TO <- stan_surv(Surv(ibi.adj, ibi.censored) ~ tve(mom.z) + tve(mom.z.2) + (1|mom.id) + (1|unit) + (1|year), 
                       data=ibi.gelada[ibi.gelada$takeover==0,], cores=2, chains=4, basehaz = 'ms')
print(summary(fit.tg.NTO, prob=c(0.025,0.975)), digits=2)

## PART 3: INFANT SURVIVAL
## RUN MODELS
fit.pu2 <- stan_surv(Surv(inf.age.out, died) ~ mom.z + mom.z.2 + (1|mom.id) + (1|year), 
                    data=surv.chacma, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.pu2, prob=c(0.025,0.975)), digits=2)

surv.gelada$inf.age.out[surv.gelada$inf.age.out==0]<-0.001
fit.tg2 <- stan_surv(Surv(inf.age.out, died) ~ mom.z + mom.z.2 + (1|mom.id) + (1|unit) + (1|year), 
                    data=surv.gelada, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.tg2, prob=c(0.025,0.975)), digits=2)

## RUN MODELS -- PARITY
surv.chacma$parity<-as.factor(surv.chacma$parity)
surv.chacma$parity<-factor(surv.chacma$parity, levels=c("Primiparous","Multiparous"))

fit.pu2.parity <- stan_surv(Surv(inf.age.out, died) ~ mom.z + parity + (1|mom.id) + (1|year), 
                     data=surv.chacma, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.pu2.parity, prob=c(0.025,0.975)), digits=2)

surv.gelada$inf.age.out[surv.gelada$inf.age.out==0]<-0.001
fit.tg2.parity <- stan_surv(Surv(inf.age.out, died) ~ mom.z + parity + (1|mom.id) + (1|unit) + (1|year), 
                     data=surv.gelada, cores=2, chains=4, basehaz = 'ms')
print(summary(fit.tg2.parity, prob=c(0.025,0.975)), digits=2)

## COMPARE CHACMA MODELS
loo1<-loo(fit.pu2)
loo2<-loo(fit.pu2.parity)
loo_compare(loo1, loo2)

## COMPARE GELADA MODELS
loo1<-loo(fit.tg2)
loo2<-loo(fit.tg2.parity)
loo_compare(loo1, loo2)

## PLOT PROJECTIONS -- CHACMA BABOONS
vals<-round(quantile(surv.chacma$mom.z, probs=c(.1,.3,.5,.7,.9)),2)
age.effect <- dplyr::select(surv.chacma, mom.id, year) %>%
  tidyr::expand(mom.z=c(vals),
                nesting(year,mom.id)) %>%
  dplyr::mutate(id = 1:n(),
                age.cats = as.factor(mom.z), mom.z.2=mom.z^2) %>%
  nest_legacy(-c(age.cats))  %>%
  dplyr::mutate(., preds = map(
    data,
    ~ posterior_survfit(
      fit.pu2, prob=0.89,
      newdata = .x,
      times = 0,
      standardise = TRUE,
      extrapolate = TRUE,
      dynamic = TRUE
    )
  ))

age.effect3 <- unnest(age.effect, preds) %>%
  dplyr::select(-data)

levels(age.effect3$age.cats)<-c("10th","30th","50th","70th","90th")

plot.pu2<-ggplot(age.effect3,
                 aes(
                   x=time, y=median,
                   color=age.cats
                 )) + labs(fill="age category") + labs(color="Maternal age\npercentile") +
  geom_line(lwd=1.25,alpha=0.8)  + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold")) +
  scale_color_manual(values=cb) +
  xlab("offspring age (years)") + ylab("proportion offspring surviving") + theme(legend.position = "none") +
  ylim(c(0.5, 1)) + scale_x_continuous(breaks=seq(0,1,0.5),labels=seq(0,1,0.5)) 
plot.pu2

## PLOT PROJECTIONS -- GELADAS
vals<-round(quantile(surv.gelada$mom.z, probs=c(.1,.3,.5,.7,.9)),2)
age.effect <- dplyr::select(surv.gelada, mom.id, year,unit) %>%
  tidyr::expand(mom.z=c(vals),
                nesting(unit, year,mom.id)) %>%
  dplyr::mutate(id = 1:n(),
                age.cats = as.factor(mom.z), mom.z.2=mom.z^2) %>%
  nest_legacy(-c(age.cats))  %>%
  dplyr::mutate(., preds = map(
    data,
    ~ posterior_survfit(
      fit.tg2, prob=0.89,
      newdata = .x,
      times = 0,
      standardise = TRUE,
      extrapolate = TRUE,
      dynamic = TRUE
    )
  ))

age.effect4 <- unnest(age.effect, preds) %>%
  dplyr::select(-data)

levels(age.effect4$age.cats)<-c("10th","30th","50th","70th","90th")

plot.tg2<-ggplot(age.effect4,
                 aes(
                   x=time, y=median,
                   color=age.cats
                 )) + labs(fill="age category") + labs(color="Maternal age\npercentile") +
  geom_line(lwd=1.25,alpha=0.8)  + theme(plot.title=element_text(face="bold", hjust=0.5), legend.title=element_text(face="bold")) +
  scale_color_manual(values=cb) +
  annotate('text', label="young", x=0.6, size=5, y=0.98) +
  annotate('text', label="old", x=0.55, size=5, y=0.81) +
  xlab("offspring age (years)") + ylab("") +
  ylim(c(0.5, 1)) + scale_x_continuous(breaks=seq(0,1,0.5),labels=seq(0,1,0.5)) 
plot.tg2
plot.pu1+plot.tg1+plot.pu2+plot.tg2 + plot_annotation(tag_levels = "a")
ggsave("Repro_aging_chacmas_v_geladas.jpg", units='cm', height=20, width=24)

## PART 4: INFANT MORTALITY CAUSES

## VISUALIZE EFFECTS -- FIGURE 4
## CHACMAS
chacma.mortality.summary1<-surv.chacma[!(surv.chacma$cause=="Censored"),] %>% 
  group_by(age.cat) %>%
  dplyr::summarise(Total.Infants=n())

sum(chacma.mortality.summary1$Total.Infants)

chacma.mortality.summary2<-surv.chacma[!(surv.chacma$cause=="Censored"),] %>% 
  group_by(age.cat, cause) %>%
  dplyr::summarise(Infants=n())

chacma.to.plot<-left_join(chacma.mortality.summary1, chacma.mortality.summary2, by="age.cat")
chacma.to.plot$Proportion<-chacma.to.plot$Infants/chacma.to.plot$Total.Infants
chacma.to.plot<-chacma.to.plot[!chacma.to.plot$cause=="Survived",]
chacma.to.plot$age.cat<-as.factor(chacma.to.plot$age.cat)
chacma.to.plot$age.cat<-factor(chacma.to.plot$age.cat, levels=c("young","mid","old"))
levels(chacma.to.plot$age.cat)<-c("Young","Middle","Old")

branded_colors <- c("#00798c","#66a182", "#440154FF", "gray")
plot1<-ggplot(data=chacma.to.plot, aes(x=age.cat, y=Proportion, fill=cause)) + geom_col(alpha=0.8) + labs(fill="probable cause") +
  scale_fill_manual(values = branded_colors) + xlab("age category") + ylab("proportion of infants") +  ggtitle("Chacma baboons") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),legend.position = "none") + ylim(0,0.5)
plot1

## GELADAS
gelada.mortality.summary1<-surv.gelada[!(surv.gelada$cause=="Censored"),] %>% 
  group_by(age.cat) %>%
  dplyr::summarise(Total.Infants=n())

sum(gelada.mortality.summary1$Total.Infants)

gelada.mortality.summary2<-surv.gelada[!(surv.gelada$cause=="Censored"),] %>% 
  group_by(age.cat, cause) %>%
  dplyr::summarise(Infants=n())

gelada.to.plot<-left_join(gelada.mortality.summary1, gelada.mortality.summary2, by="age.cat")
gelada.to.plot$Proportion<-gelada.to.plot$Infants/gelada.to.plot$Total.Infants

gelada.to.plot<-gelada.to.plot[!gelada.to.plot$cause=="Survived",]
gelada.to.plot$age.cat<-as.factor(gelada.to.plot$age.cat)
gelada.to.plot$age.cat<-factor(gelada.to.plot$age.cat, levels=c("young","mid","old"))
levels(gelada.to.plot$age.cat)<-c("Young","Middle","Old")

plot2<-ggplot(data=gelada.to.plot, aes(x=age.cat, y=Proportion,fill=cause)) + geom_col(alpha=0.8) + labs(fill="probable cause") +
  scale_fill_manual(values = branded_colors) + xlab("age category") + ylab("") +  ggtitle("Geladas") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) + ylim(0,0.5)
plot1+plot2 +plot_annotation(tag_levels = "a")
ggsave(file="Infant mortality_causes.jpg", units="cm", width=21, height=10, dpi=300)

## MAKE MODELS
## SET INFANTICIDES AND DEATHS
surv.chacma$Infanticide<-0
surv.chacma$Infanticide[surv.chacma$cause=="Infanticide"]<-1
surv.gelada$Infanticide<-0
surv.gelada$Infanticide[surv.gelada$cause=="Infanticide"]<-1

surv.chacma$MatDeath<-0
surv.chacma$MatDeath[surv.chacma$cause=="Maternal death"]<-1
surv.gelada$MatDeath<-0
surv.gelada$MatDeath[surv.gelada$cause=="Maternal death"]<-1

surv.chacma$Injured<-0
surv.chacma$Injured[surv.chacma$cause=="Injury or illness"]<-1
surv.gelada$Injured<-0
surv.gelada$Injured[surv.gelada$cause=="Injury or illness"]<-1

surv.chacma$Unknown<-0
surv.chacma$Unknown[surv.chacma$cause=="Unknown"]<-1
surv.gelada$Unknown<-0
surv.gelada$Unknown[surv.gelada$cause=="Unknown"]<-1

## REMOVE CENSORED CASES
chacma.causes<-surv.chacma[!surv.chacma$cause=="Censored",]
gelada.causes<-surv.gelada[!surv.gelada$cause=="Censored",]

## INFANTICIDE
infanticide.chacma<-brm(data=chacma.causes, Infanticide ~ mom.z + mom.z.2 + (1|year), 
            prior=c(
              prior(normal(0, 5), "Intercept"),
              prior(normal(0, 1), "b")),
            family="bernoulli", chains=4, core=2)
summary(infanticide.chacma)

infanticide.gelada<-brm(data=gelada.causes, Infanticide ~ mom.z + mom.z.2 + (1|year), 
                        prior=c(
                          prior(normal(0, 5), "Intercept"),
                          prior(normal(0, 1), "b")),
                        family="bernoulli", chains=4, core=2)
summary(infanticide.gelada)

## MATERNAL DEATH
mat.chacma<-brm(data=chacma.causes, MatDeath ~ mom.z + mom.z.2 + (1|year), 
            prior=c(
              prior(normal(0, 5), "Intercept"),
              prior(normal(0, 1), "b")),
            family="bernoulli", chains=4, core=2)
summary(mat.chacma)

mat.gelada<-brm(data=gelada.causes, MatDeath ~ mom.z + mom.z.2 + (1|year), 
                        prior=c(
                          prior(normal(0, 5), "Intercept"),
                          prior(normal(0, 1), "b")),
                        family="bernoulli", chains=4, core=2)
summary(mat.gelada)

## INJURY OR ILLNESS
inj.chacma<-brm(data=chacma.causes, Injured ~ mom.z + mom.z.2 + (1|year), 
                prior=c(
                  prior(normal(0, 5), "Intercept"),
                  prior(normal(0, 1), "b")),
                family="bernoulli", chains=4, core=2)
summary(inj.chacma)

inj.gelada<-brm(data=gelada.causes, Injured ~ mom.z + mom.z.2 + (1|year), 
                prior=c(
                  prior(normal(0, 5), "Intercept"),
                  prior(normal(0, 1), "b")),
                family="bernoulli", chains=4, core=2)
summary(inj.gelada)

## UNKNOWN
unk.chacma<-brm(data=chacma.causes, Unknown ~ mom.z + mom.z.2 + (1|year), 
                prior=c(
                  prior(normal(0, 5), "Intercept"),
                  prior(normal(0, 1), "b")),
                family="bernoulli", chains=4, core=2)
summary(unk.chacma)

unk.gelada<-brm(data=gelada.causes, Unknown ~ mom.z + mom.z.2 + (1|year), 
                prior=c(
                  prior(normal(0, 5), "Intercept"),
                  prior(normal(0, 1), "b")),
                family="bernoulli", chains=4, core=2)
summary(unk.gelada)

