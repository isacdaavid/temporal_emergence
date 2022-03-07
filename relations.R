library(ggplot2)
theme_set(theme_gray(base_size = 22))

substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
}

df <- read.csv('./relations.csv', header = TRUE, sep = ",")
df$neurons_state <- paste(df$neurons, df$state, sep='_')

## big Phi, connected vs disconnected
dftemp  <- data.frame(
    'neurons'=sapply(unique(df$neurons_state), function(x) {unique(df[df$neurons_state == x, 'neurons'])}),
    'state'=sapply(unique(df$neurons_state), function(x) {unique(df[df$neurons_state == x, 'state'])}),
    'neurons_state'=unique(df$neurons_state),
    'Φ'=sapply(unique(df$neurons_state), function(x) {unique(df[df$neurons_state == x, 'big_Phi'])}),
    'Connectivity'=sapply(unique(df$neurons_state), function(x) {unique(df[df$neurons_state == x, 'connected'])})
)

dftemp$Connectivity  <- sapply(dftemp$Connectivity,
                               function(x) {if (x == 'True') 'Connected' else 'Disconnected'})
## dftemp$Φ  <- sapply(dftemp$Φ, function (phi) {if (phi <= 0) 1e-14 else phi})

ggplot(dftemp, aes(x=Connectivity,
                   y=Φ,
                   fill=Connectivity,
                   color=Connectivity)) +
    geom_violin(alpha=.5) +
    geom_jitter(width=.1) +
    labs(title='Big Φ per neuron pair',
         y='log(Φ)') +
    guides(fill='none', color='none') +
    scale_y_log10()

t.test(
    sapply(dftemp[dftemp$Connectivity == 'Connected', 'Φ'], log),
    sapply(dftemp[dftemp$Connectivity == 'Disconnected', 'Φ'], log),
    alternative='greater',
    var.equal = FALSE
)

## amount of relations, connected vs disconnected
y <- sapply(unique(df$neurons), function(pair) {nrow(df[df$neurons == pair,])})
x <- sapply(
    sapply(unique(df$neurons), function(pair) {unique(df[df$neurons == pair, 'connected'])}),
    function(x) {if (x == 'True') 'Connected' else 'Disconnected'}
    )

dftemp2 <- data.frame('Amount_of_relations' = y, 'Connectivity' = x)

ggplot(dftemp2, aes(y=Amount_of_relations,
                   x=Connectivity,
                   fill=Connectivity,
                   color=Connectivity)) +
    geom_violin(alpha=.5) +
    geom_jitter(width=.1) +
    labs(title='Amount of relations per neuron pair',
         y='log(Amount of relations)') +
    guides(fill='none', color='none') +
    scale_y_log10()

t.test(
    sapply(dftemp2[dftemp2$Connectivity == 'Connected', 'Amount_of_relations'], log),
    sapply(dftemp2[dftemp2$Connectivity == 'Disconnected', 'Amount_of_relations'], log),
    alternative='greater',
    var.equal = FALSE
    )

## amount of relations, connected vs disconnected (by state)
dftemp3 <- df[, c('neurons', 'state', 'big_Phi', 'connected')]
dftemp3$neurons_state  <- paste(dftemp3$neurons, dftemp3$state, sep= '_')
y <- sapply(unique(dftemp3$neurons_state), function(x) {nrow( dftemp3[dftemp3$neurons_state == x, ] )})
x  <- sapply(
    unique(dftemp3$neurons_state),
    function(x) {unique(dftemp3[dftemp3$neurons_state == x, 'connected'])}
)
x  <- sapply(x, function(x) {if (x == 'True') 'Connected' else 'Disconnected'})
x2 <- sapply(names(x), function(cosa) {substrRight(cosa, 6)})
dftemp3.2 <- data.frame('Amount_of_relations' = y,
                        'State' = x2,
                        'Connectivity' = x)

ggplot(dftemp3.2, aes(y=Amount_of_relations,
                   x=State,
                   fill=Connectivity,
                   color=Connectivity)) +
    geom_jitter(width=.2) +
    geom_violin(alpha=.5) +
    labs(title='Amount of relations per neuron pair',
         y='log(Amount of relations)') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_log10()

dftemp3.2$log_amount  <- sapply(dftemp3.2$Amount_of_relations, log)
model <- aov(log_amount ~ State * Connectivity, data = dftemp3.2)
summary(model)
hist(model$residuals, main = "Histogram of Residuals", xlab = "Residuals", col = "steelblue")
TukeyHSD(model, conf.level=.95)

## correlation between phi and # or relations
dftemp3.2$neurons_state <- row.names(dftemp3.2)
dftemp4 <- data.frame(
    Φ=sapply(dftemp$neurons_state, function(x) {
        phi <- dftemp[dftemp$neurons_state == x, 'Φ']
        if (phi < 0 || phi < 1e-14) 0 else phi
    }),
    Amount_of_relations=sapply(
        dftemp$neurons_state,
        function(x) {dftemp3.2[dftemp3.2$neurons_state == x, 'Amount_of_relations']}
    ),
    Connectivity=sapply(dftemp$neurons_state, function(x) {dftemp[dftemp$neurons_state == x, 'Connectivity']})
)

kk  <- data.frame(x=sapply(dftemp4$Φ, function(x) {
    x2  <- if (x <= 0) 1e-14 else x
    log(x2)
}),
y=sapply(dftemp4$Amount_of_relations, function(x) {log(x)}))

model  <- lm(y ~ x, kk)

ggplot(dftemp4, aes(x=Φ, y=Amount_of_relations, color=Connectivity)) +
    geom_jitter() +
    geom_smooth(method='lm', formula=y~x) +
    scale_y_log10() +
    scale_x_log10() +
    labs(title="Poisson regression",
         x='log(Φ)',
         y='log(Amount of relations)')
