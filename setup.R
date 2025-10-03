install.packages("renv")
renv::init()

renv::settings$snapshot.type("all")
renv::settings$package.dependency.fields(c("Imports","Depends","LinkingTo","Suggests"))

sink("sessionInfo.txt"); sessionInfo(); sink()
