
pain = c(4, 5, 4, 3, 2, 4, 3, 4, 4, 6, 8, 4, 5, 4, 6, 5, 8, 6, 6, 7, 6, 6, 7, 5, 6, 5, 5)
drug = c(rep("A",9), rep("B",9), rep("C",9))
migraine = data.frame(pain,drug)
migraine


plot(pain ~ drug, data=migraine)


results = aov(pain ~ drug, data=migraine)
summary(results)



pairwise.t.test(x = pain, g = drug, p.adjust="bonferroni") # g is the grouping vector


# Another multiple comparisons procedure is Tukey‟s method (a.k.a. Tukey's Honest Significance Test). The function TukeyHSD() creates a set of confidence intervals on the differences between means with the specified family-wise probability of coverage. The general form is
# TukeyHSD(x, conf.level = 0.95)
# Here x is a fitted model object (e.g., an aov fit) and conf.level is the confidence level.
# TukeyHSD(results, conf.level = 0.95)

TukeyHSD(results, conf.level = 0.95)



