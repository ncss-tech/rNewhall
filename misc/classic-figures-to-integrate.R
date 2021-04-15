



p1 = ggplot(soil.profile.matrix.plot, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  labs(x="", y="",
       title=paste("Month ", m, sep = "")) +
  scale_y_reverse(breaks=c(1:8),labels = c("Surface", "Moisture Control", "Moisture Control","","","","","Bottom"))+
  scale_x_continuous(breaks=c(1:8),labels = c("PWP", "", "","","","","","AWC")) +
  scale_fill_gradient(limits = c(0,1)) +#, colours=c("#132B43","#56B1F7" ))+
  theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=10),
                     plot.title=element_text(size=14),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  geom_hline(yintercept=.5) + geom_hline(yintercept=1.5) + geom_hline(yintercept=2.5) + geom_hline(yintercept=3.5) + geom_hline(yintercept=4.5) +
  geom_hline(yintercept=5.5) + geom_hline(yintercept=6.5) + geom_hline(yintercept=7.5) + geom_hline(yintercept=8.5) + geom_vline(xintercept=.5) +
  geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) + geom_vline(xintercept=4.5) + geom_vline(xintercept=5.5) + geom_vline(xintercept=6.5) +
  geom_vline(xintercept=7.5) + geom_vline(xintercept=8.5) +labs(fill = "% Filled")


p2 = ggplot(soil.profile.matrix.plot, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  labs(x="", y="",
       title=paste("Month ", m, sep = "")) +
  scale_y_reverse(breaks=c(1:8),labels = c("Surface", "Moisture Control", "Moisture Control","","","","","Bottom"))+
  scale_x_continuous(breaks=c(1:8),labels = c("PWP", "", "","","","","","AWC")) +
  scale_fill_gradient(limits = c(0,1))+#, colours=c("#132B43","#56B1F7" ))+
  theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=10),
                     plot.title=element_text(size=14),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  geom_hline(yintercept=.5) + geom_hline(yintercept=1.5) + geom_hline(yintercept=2.5) + geom_hline(yintercept=3.5) + geom_hline(yintercept=4.5) +
  geom_hline(yintercept=5.5) + geom_hline(yintercept=6.5) + geom_hline(yintercept=7.5) + geom_hline(yintercept=8.5) + geom_vline(xintercept=.5) +
  geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) + geom_vline(xintercept=4.5) + geom_vline(xintercept=5.5) + geom_vline(xintercept=6.5) +
  geom_vline(xintercept=7.5) + geom_vline(xintercept=8.5)+labs(fill = "% Filled")


p3 = ggplot(soil.profile.matrix.plot, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  labs(x="", y="",
       title=paste("Month ", m, sep = "")) +
  scale_y_reverse(breaks=c(1:8),labels = c("Surface", "Moisture Control", "Moisture Control","","","","","Bottom"))+
  scale_x_continuous(breaks=c(1:8),labels = c("PWP", "", "","","","","","AWC")) +
  scale_fill_gradient(limits = c(0,1))+#, colours=c("#132B43","#56B1F7" ))+
  theme_bw() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=10),
                     plot.title=element_text(size=14),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  geom_hline(yintercept=.5) + geom_hline(yintercept=1.5) + geom_hline(yintercept=2.5) + geom_hline(yintercept=3.5) + geom_hline(yintercept=4.5) +
  geom_hline(yintercept=5.5) + geom_hline(yintercept=6.5) + geom_hline(yintercept=7.5) + geom_hline(yintercept=8.5) + geom_vline(xintercept=.5) +
  geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) + geom_vline(xintercept=4.5) + geom_vline(xintercept=5.5) + geom_vline(xintercept=6.5) +
  geom_vline(xintercept=7.5) + geom_vline(xintercept=8.5) +labs(fill = "% Filled")
