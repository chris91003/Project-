#Iterate through data frame and take each column and convert to Vector 
#iterate through columns and maybe see if there are GO terms that are in multiple clusters?
#iterate using other tissues????


df <- read.csv("kendall-liver-gene-clusters-high-variance.csv")
df

View(df)

#Create empty vector that will contain gene cluster columns
column_vector <- vector()

#for loop to iterate through df and extract columns and convert to vector 

for (i in 1:ncol(df)) {
  vector <- c(column_vector, df)
  print(vector)
  }

