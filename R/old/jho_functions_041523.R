### Functions

### Updated 041524
my_heatmaply_cor <- function(comm){
    r <- cor(comm)
    ## not mine, from heatmaply docs
    ## We use this function to calculate a matrix of p-values from correlation tests
    ## https://stackoverflow.com/a/13112337/4747043
    cor.test.p <- function(x){
        FUN <- function(x, y) cor.test(x, y)[["p.value"]]
        z <- outer(
            colnames(x), 
            colnames(x), 
            Vectorize(function(i,j) FUN(x[,i], x[,j]))
        )
        dimnames(z) <- list(colnames(x), colnames(x))
        z
    }
    p <- cor.test.p(comm)
    p_note <- (p < 0.05) %>% as.numeric() %>% matrix(nrow = ncol(comm), ncol = ncol(comm))
    
    heatmaply_cor(
        r*p_note,
        node_type = "scatter",
        point_size_mat = -log(p) * p_note, 
        point_size_name = "-log(p-value)",
        label_names = c("x", "y", "Correlation")
    ) %>% return()
}

my_ggsave <- function(plot, 
                      fig_num,
                      filename = NA,
                      device = 'png',
                      path = './Figs/', 
                      dpi = 160, 
                      units = 'mm', 
                      width = 120, 
                      height = 100, 
                      bg = 'white'){
    library(ggpubr)
    if (is.na(filename)) {
        filename = paste('figure_', as.character(fig_num), '.', device, sep = '')
    }
    ggsave(
        filename,
        plot,
        device = device,
        path = path,
        dpi = dpi,
        units = units,
        width = width,
        height = height,
        bg = bg
    )
}

tiles <- function(mtrx, mid = 0.05){
    if(is.matrix(mtrx)){
        mtrx <- melt(mtrx)
    }
    ggplot(mtrx, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(na.rm = FALSE) +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_fill_gradient2(low = "red", high = "green2", midpoint = mid, mid = "gray95", na.value = "white") +
        geom_text(aes(label = round(value, digits = 1)), size = 2)
}

p3d <- function(df, 
                axis1, axis2, axis3,
                breaks = 10, 
                phi = 15, 
                shape = NA, 
                color = c(), 
                midpoint = NA,
                display = FALSE,
                save = FALSE,
                filename = "3d_plot",
                extension = ".png",
                path = str_c("Plots/", today(), "/"),
                vector = NULL){
    p3d_list <- list()
    for(i in 1:(360/breaks)){
        theta <- breaks*i
        if(is.na(midpoint)){midpoint <- median(color)}
        p <- ggplot(data = df, aes(x = axis1, y = axis2, z = axis3,
                                   shape = shape, color = color)) +
            axes_3D(phi = phi, theta = theta) +
            stat_3D(phi = phi, theta = theta) +
            labs_3D(phi = phi, theta = theta) +
            scale_color_gradient2(low = "black",
                                  midpoint = midpoint, 
                                  mid = "cyan2", 
                                  high = "green") +
            labs(color = "Gradient", shape = "Factor",
                 x = element_blank(), y = element_blank()) %>%
            suppressWarnings()
        # can accept one ordination vector
        if(length(vector)){
            p <- p + stat_3D(data = vector, phi = 10, theta = theta, geom = "line",
                             aes(x = vector[,1], y = vector[,2], z = vector[,3],
                                 color = NULL, shape = NULL, ), arrow = arrow())
        }
        if(save){
            ggsave(str_c(path, filename, "/", filename,
                         str_pad(as.character(i), width = 3, pad = "0"),
                         extension),
                   width = 1600, height = 1000, units = "px")
        }
        p3d_list[[i]] <- p
    }
    if(save){
        ## list file names and read in
        imgs <- list.files(str_c(path, "/", filename, "/"), full.names = TRUE)
        img_list <- lapply(imgs, image_read) 
        img_joined <- image_join(img_list) # join the images together
        img_animated <- image_animate(img_joined, fps = .5) # animate
        if(display){print(img_animated)} # view animated image
        image_write(image = img_animated,
                    path = str_c(path, filename, extension = ".gif")) # save to disk
        return(img_animated)
    }
    
    return(p3d_list)
}

# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    par(usr = c(0, 1, 0, 1))
    r2 <- (cor(x, y))^2
    txt <- format(c(r2, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r2^.25) ###
}

ses.plot <- function(pbar, xlim = c(0, 160)){ # takes permustats() %>% summary() object
    pbar <- pbar$z
    pbar <- pbar %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        rename(Factor = rowname, Stat = ".")
    plot <- ggplot(pbar, aes(y = reorder(Factor, -Stat), x = Stat)) + 
        geom_bar(stat = 'identity', position = 'stack', fill = 'cyan4') +
        xlab("SES after marginal PerMANOVA") +
        coord_cartesian(xlim = xlim) +
        ggtitle("Residualized effects on community") +
        theme_few(base_size = 18) +
        theme(axis.title.y = element_blank())
    return(plot)
}

p2d <- function(df, shape = '', color = '', display = "sites", midpoint = 100000){
    if(display == 'sites'){
        obj <- metaMDS(df)$points %>%
            as.data.frame()
    }else if(display == 'species'){
        obj <- metaMDS(df)$species %>%
            as.data.frame()
    }else{
        return(errorCondition('invalid display'))
    }
    plot <- ggplot(obj, aes(x = MDS1, y = MDS2, shape = shape, color = color, fill = color)) + 
        scale_fill_gradient2(low = "black",
                             midpoint = midpoint, 
                             mid = "cyan2", 
                             high = "green") +
        scale_color_gradient2(low = "black",
                              midpoint = midpoint, 
                              mid = "cyan2", 
                              high = "green") +
        scale_shape_manual(values=c(21, 22, 23, 24, 25, 8, 9)) +
        geom_point(show.legend = T) +
        stat_ellipse(aes(group = shape), type = 'euclid', level = .5) +
        ggtitle("nMDS of grouped observations") +
        theme_few(base_size = 18) + 
        labs(shape = "Shape", color = "Color", fill = "Color")
    
    return(plot)
}

mds_jho <- function(c, m, k = 2, trymax = 20){
    mds <- metaMDS(c, k = k, trymax = trymax)
    mds$points %>%
        as.data.frame() %>%
        rownames_to_column('id') %>%
        full_join(rownames_to_column(m, 'id'), by = 'id') %>%
        column_to_rownames('id') %>%
        return()
}

pcoa_jho <- function(c, m){
    comm_pcoa <- vegdist(c, method = 'bray') %>%
        pcoa()
    vectors <- comm_pcoa$vectors %>%
        as.data.frame() %>%
        rownames_to_column('id') %>%
        full_join(rownames_to_column(m, 'id'), by = 'id') %>%
        column_to_rownames('id') %>%
        return()
}

betadisper_plot <- function(bd, t){
    vect <- as.data.frame(bd$vectors)
    vect$t <- as.character(t)
    dist <- as.data.frame(bd$group.distances) %>%
        rownames_to_column(var = 't')
    centroids <- as.data.frame(bd$centroids) %>% 
        rownames_to_column(var = 't') %>%
        select(PCoA1,PCoA2,t) %>%
        rename(c1 = PCoA1, c2 = PCoA2)
    vectj <- left_join(vect, centroids)
    vectj <- left_join(vectj, dist)
    p <- ggplot(vectj, aes(x = PCoA1, y = PCoA2, shape = t, fill = t, group = t, color = t)) +
        scale_shape_manual(values=c(21, 22, 23, 24, 25, 8, 9)) +
        geom_point(show.legend = T) +
        stat_ellipse(level = x, show.legend = F) +
        geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
        geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
        geom_segment(aes(x=c1, y=c2, xend=PCoA1, yend=PCoA2), 
                     size=1, show.legend=FALSE) +
        # geom_point(data=sites.long2, 
        #        aes(x=axis1, y=axis2, colour=Management, shape=Management), 
        #        size=5) +
        
        ggtitle("PCoA of grouped betadispersions") +
        theme_few(base_size = 18)
    
    return(p)
}

nmds_follow_score <- function(df, df2){ # accepts dfs with MDS1/MDS2 vars
    
    df$MDS1c <- df$MDS1
    df$MDS2c <- df$MDS2
    df$MDS1e <- df2$MDS1
    df$MDS2e <- df2$MDS2
    
    df$follow_score <- sqrt((df$MDS2c - df$MDS1c)^2 + (df$MDS2e - df$MDS1e)^2)
    
    return(df$follow_score) # returns vector with line lengths
}

# modified from FD package
jfunctcomp <- function (x, a, CWM.type = c("dom", "all"), bin.num = NULL)
{
    # From package 'FD'
    # editing for my use
    # removed 'dom' dominant cwm functionality
    if(CWM.type == 'dom'){
        print("I've removed this functionality bc it was throwing errors and I don't need it.")
        return(0)
    }
    if (!is.matrix(x) & !is.data.frame(x))
        stop("'x' must be a matrix or a data frame.", "\n")
    else x <- data.frame(x)
    if (!is.matrix(a))
        stop("'a' must be a matrix.", "\n")
    if (is.null(row.names(x)))
        stop("'x' must have row names.", "\n")
    else x.n <- row.names(x)
    if (is.null(colnames(a)))
        stop("'a' must have column names.", "\n")
    else a.n <- colnames(a)
    s.x <- dim(x)[1]
    s.a <- dim(a)[2]
    if (s.x != s.a)
        stop("Different number of species in 'x' and 'a'.",
            "\n")
    if (any(x.n != a.n))
        stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).",
            "\n")
    com <- dim(a)[1]
    t <- dim(x)[2]
    com.names <- row.names(a)
    sp.names <- row.names(x)
    tr.names <- names(x)
    a[which(is.na(a))] <- 0
    CWM.type <- match.arg(CWM.type)
    is.bin <- function(k) all(k[!is.na(k)] %in% c(0, 1))
    bin.var <- rep(NA, t)
    names(bin.var) <- tr.names
    for (i in 1:t) bin.var[i] <- is.bin(x[, i])
    if (!all(bin.var[bin.num]))
        stop("'bin.num' points to non-binary variables.\n")
    bin.var[bin.num] <- FALSE
    type <- sapply(x, data.class)
    type[type %in% c("numeric", "integer")] <- "C"
    type[type == "ordered"] <- "O"
    type[type == "factor"] <- "N"
    type[bin.var] <- "B"
    sum.a <- apply(a, 1, sum)
    a <- a/sum.a
    a <- t(a)
    a <- data.frame(a)
    temp <- list()
    for (i in 1:t) {
        if (type[i] == "C") {
            vec <- numeric(com)
            for (j in 1:com) vec[j] <- weighted.mean(x[, i],
                a[, j], na.rm = T)
            temp[[i]] <- matrix(vec, com, 1, dimnames = list(com.names,
                tr.names[i]))
        }
        else {
            x[, i] <- as.factor(x[, i])
            fac <- data.frame()
            which.dom <- rep(NA, com)
            for (k in 1:com) {
                temp2 <- tapply(a[, k], x[, i], sum)
                fac <- rbind(fac, temp2)
                #which.dom[k] <- sample(levels(x[, i])[which(fac[k,] == max(fac[k,]))], size = 1)
            }
            colnames(fac) <- paste(tr.names[i], "_", levels(x[,
                i]), sep = "")
            rownames(fac) <- com.names
            #which.dom <- data.frame(which.dom)
            #colnames(which.dom) <- tr.names[i]
            #rownames(which.dom) <- com.names
            #if (CWM.type == "dom")
            #    temp[[i]] <- which.dom
            if (CWM.type == "all")
                temp[[i]] <- fac
        }
    }
    temp <- data.frame(temp)
    return(temp)
}

