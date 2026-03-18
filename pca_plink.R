setwd('~/Desktop/AF')
data = read.table('AF_pca.eigenvec')
popinfo = read.table('sample_popinfo.tsv',header = T)
colnames(data) <- c('fid', 'iid', 
                    'pc1', 'pc2', 'pc3','pc4', 'pc5',
                    'pc6','pc7', 'pc8', 'pc9','pc10')
data$fid = popinfo$Region

library(plotly)
plot_ly(data, 
        x = ~pc1, 
        y = ~pc2, 
        z = ~pc3, 
        color = ~fid,     
        colors = "Set1",  
        type = 'scatter3d', 
        mode = 'markers',
        marker = list(size = 4, opacity = 0.8),
        text = ~iid)

