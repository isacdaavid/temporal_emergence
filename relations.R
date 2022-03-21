set.seed(1234)
library(ggplot2)
theme_set(theme_gray(base_size = 22))

FIG_PATH <- './figures/relations/relations.R'

substrRight <- function(x, n) {
    substr(x, nchar(x)-n+1, nchar(x))
}

df <- read.csv('./relations.csv', header = TRUE, sep = ",")
df$neurons_state <- paste(df$neurons, df$state, sep='_')

## discard rows with negative Big Phi
## df  <- df[df$big_Phi >= 0, ]

######################################################################
## big Phi, connected vs disconnected
######################################################################

dftemp  <- data.frame(
    'neurons'=sapply(unique(df$neurons_state),
                     function(x) {unique(df[df$neurons_state == x, 'neurons'])}),
    'state'=sapply(unique(df$neurons_state),
                   function(x) {unique(df[df$neurons_state == x, 'state'])}),
    'neurons_state'=unique(df$neurons_state),
    'Φ'=sapply(unique(df$neurons_state),
               function(x) {unique(df[df$neurons_state == x, 'big_Phi'])}),
    'Connectivity'=sapply(unique(df$neurons_state),
                          function(x) {
                              val <- unique(df[df$neurons_state == x, 'connected'])
                              if (val == 'True') 'Connected' else 'Disconnected'
                          })
)

wilcox.test(Φ ~ Connectivity, data=dftemp, alternative='greater')

svg(paste0(FIG_PATH, '/Φ_conn_vs_disconn.svg'))
ggplot(dftemp, aes(x=Connectivity,
                   y=Φ,
                   fill=Connectivity,
                   color=Connectivity)) +
    geom_violin(alpha=.2) +
    geom_jitter(aes(shape=state), width=.4, alpha=.5) +
    labs(title='Big Φ per system',
         y='Φ') +
    guides(fill='none', color='none')
dev.off()

######################################################################
## amount of relations, connected vs disconnected
######################################################################

dftemp$Amount_of_relations <- sapply(
    unique(df$neurons_state),
    function(x) {nrow(df[df$neurons_state == x, ])}
)

## normalize by number of possible relations
## 2^(2*n) - 2*n - 1 : 2^m - m - 1
dftemp$Amount_of_possible_relations  <- sapply(
    unique(df$neurons_state),
    function(x) {
        mechs <- unique(df[df$neurons_state == x, 'system_mechs'])
        2^(2*mechs) - 2*mechs - 1
    }
)

dftemp$Amount_of_relations_norm <- dftemp$Amount_of_relations / dftemp$Amount_of_possible_relations

wilcox.test(Amount_of_relations_norm ~ Connectivity, data=dftemp, alternative='greater')

svg(paste0(FIG_PATH, '/amount_relations_norm_conn_vs_disconn.svg'))
ggplot(dftemp, aes(y=Amount_of_relations_norm,
                   x=Connectivity,
                   fill=Connectivity,
                   color=Connectivity)) +
    geom_violin(alpha=.2) +
    geom_jitter(mapping=aes(shape=state), width=.4, alpha=.5) +
    labs(title='Proportion of relations per system',
         y='Completeness') +
    guides(fill='none', color='none')
dev.off()

wilcox.test(Amount_of_relations ~ Connectivity, data=dftemp, alternative='greater')

svg(paste0(FIG_PATH, '/amount_relations_conn_vs_disconn.svg'))
plot <- ggplot(dftemp, aes(y=Amount_of_relations,
                           x=Connectivity,
                           fill=Connectivity,
                           color=Connectivity))

for (i in unique(dftemp$Amount_of_possible_relations)) {
    plot <- plot + geom_abline(slope=0, intercept=i, alpha=.5)
}

plot <- plot + geom_violin(alpha=.2) +
    geom_jitter(aes(shape=state), width=.4, alpha=.5) +
    labs(title='Amount of relations per system', y='Amount of relations') +
    guides(fill='none', color='none')

plot
dev.off()

######################################################################
## amount of relations, connected vs disconnected (by state)
######################################################################

svg(paste0(FIG_PATH, '/amount_relations_anova.svg'))
ggplot(dftemp, aes(y=Amount_of_relations,
                            x=state,
                            fill=Connectivity,
                            color=Connectivity)) +
    geom_jitter(aes(shape=state), width=.2, alpha=.5) +
    geom_violin(alpha=.5) +
    labs(title='Amount of relations per system',
         y='Amount of relations') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(shape='none')
dev.off()

svg(paste0(FIG_PATH, '/amount_relations_anova_norm.svg'))
ggplot(dftemp, aes(y=Amount_of_relations_norm,
                   x=state,
                   fill=Connectivity,
                   color=Connectivity)) +
    geom_jitter(aes(shape=state), width=.2, alpha=.5) +
    geom_violin(alpha=.5) +
    labs(title='Proportion of relations per system',
         y='Completeness') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(shape='none')
dev.off()

dftemp$log_amount  <- sapply(dftemp$Amount_of_relations, log)
model <- aov(log_amount ~ state * Connectivity, data = dftemp)
summary(model)
hist(model$residuals, main = "Histogram of Residuals", xlab = "Residuals", col = "steelblue")
TukeyHSD(model, conf.level=.95)

######################################################################
## correlation between phi and # or relations
######################################################################

## all data
dftemp$log_Amount_of_relations_norm  <- sapply(dftemp$Amount_of_relations_norm, log)
model  <- lm(log_Amount_of_relations_norm ~ Φ, dftemp)
summary(model)

## connected
dftemp$log_Amount_of_relations_norm  <- sapply(dftemp$Amount_of_relations_norm, log)
model  <- lm(log_Amount_of_relations_norm ~ Φ, dftemp[dftemp$Connectivity == 'Connected',])
summary(model)

## disconnected
dftemp$log_Amount_of_relations_norm  <- sapply(dftemp$Amount_of_relations_norm, log)
model  <- lm(log_Amount_of_relations_norm ~ Φ, dftemp[dftemp$Connectivity == 'Disconnected',])
summary(model)

svg(paste0(FIG_PATH, '/corr_Φ_amount_relations.svg'), width=13, height=8)
ggplot(dftemp, aes(x=Φ, y=log_Amount_of_relations_norm, color=Connectivity)) +
    geom_jitter(aes(shape=state), alpha=.75) +
    geom_smooth(method='lm', formula=y~x) +
    ## scale_y_log10() +
    ## scale_x_log10() +
    labs(title="Poisson regression",
         x='Φ',
         y='log(Completeness)')
dev.off()
