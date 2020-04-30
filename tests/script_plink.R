library(opalr)
o <- opal.login(url="https://opal-test.obiba.org")
opal.assign.resource(o, "client", "test.obesity_plink_remote")
# verify it is a SSH resource
opal.execute(o, "class(client)")
# verify allowed commands
opal.execute(o, "client$getAllowedCommands()")
# list data files
opal.execute(o, "client$exec('ls')")
# use connection ID to prefix result files
opal.execute(o, "assign('prefix', basename(tempdir()))")
# execute plink1
opal.execute(o, "client$exec('plink1', c('--bfile', 'obesity', '--assoc', '--pheno', 'obesity.phe', '--pheno-name', 'obese', '--out', paste0('/tmp/', prefix, '-out'), '--noweb'))")
# list result files
opal.execute(o, "client$exec('ls', '/tmp')")
# download result file in R server
opal.execute(o, "client$downloadFile(paste0('/tmp/', prefix, '-out.assoc'))")
# verify result file was downloaded in R server
opal.execute(o, "list.files()")
# remote server clean up
opal.execute(o, "client$exec('rm', paste0('/tmp/', prefix, '-out.*'))")
opal.execute(o, "client$close()")
opal.logout(o)

#### Debug

library(resourcer)
brge_plink <- resourcer::newResource(url="ssh://plink-test.obiba.org:2222/home/master/brge?exec=ls,plink,plink1", identity = "master", secret = "master")
client <- resourcer::newResourceClient(brge_plink)
class(client)
client$getAllowedCommands()
client$exec("ls", "-la")
tempDir <- client$tempDir()
tempDir
client$exec("ls", tempDir)
client$exec('plink1', c('--bfile', 'brge', '--assoc', '--out', paste0(tempDir, '/out1'), '--noweb'))
client$exec('plink', c('--bfile', 'brge', '--assoc', '--out', paste0(tempDir, '/out')))
client$exec("ls", tempDir)
client$removeTempDir()
client$close()







#### Create new (DSlite)
library(resourcer)
library(DSLite)
library(dsBaseClient)
library(dsBase)


dslite.server <- newDSLiteServer(resources = list(brge_plink = resourcer::newResource(url="ssh://plink-test.obiba.org:2222/home/master/brge?exec=ls,plink,plink1", identity = "master", secret = "master")))

# build login details
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "dslite.server", 
               resource = "brge_plink", driver = "DSLiteDriver")
logindata <- builder$build()# login and assign resources
o <- datashield.login(logins = logindata, assign = TRUE, symbol = "client")


ds.class("client")

dslite.server$aggregateMethod("plinkDS", plinkDS)


plink.arguments <- c("--bfile brge --freq")
ans <- ds.PLINK("client", plink.arguments)

head(ans$study1$results)

plink.arguments <- c("--bfile brge --logistic --covar brge.phe --covar-name gender,age")
ans <- ds.PLINK("client", plink.arguments)

head(ans$study1$results)


datashield.logout(o)




args <- tolist(plink.arguments)

paste0("plinkDS(", client, ", ", "'", gsub(" ''=", " ", paste(names(args), args, collapse = "', '", sep = "'='")), "')")

plinkDS("client", bfile="brge", "freq")



plink.arguments <- c("--bfile brge --freq")
gg <- unlist(strsplit(plink.arguments, " "))
gg

plinkDS("client", bfile="brge", "freq")
nn<- names(dots)
plink.command <- NULL
for (i in 1:length(nn))
{
  if (nn[i]!=""){
    arg.i <- c(paste0("--", nn[i]), dots[[i]])
    plink.command <- c(plink.command, arg.i)
  } else{
    arg.i <- paste0("--", dots[[i]])
    plink.command <- c(plink.command, arg.i)
  }
}
plink.command

ans <- ds.PLINK("client", plink.command)
lapply(ans, names)
lapply(ans$study1$results, names)
ans$study1$results

plink.command <- c("--bfile brge --freq --hardy")
ans <- ds.PLINK("client", plink.command)
lapply(ans, names)
lapply(ans$study1$results, names)





plink.command <- unlist(strsplit(plink.command, " "))
plink.command <- c(plink.command, "--noweb", "--out")
plink.command
file.out <- paste0('/tmp/', basename(tempdir()), '-out')
command <- gsub(" ", "', '", paste0("'", plink.command, "'"))
command <- paste(command, c('--out', baseout, '--noweb'))

client$exec('plink1', command)
client$downloadFile(paste0('/tmp/', prefix, '-out.assoc'))
  


#### DataSHIELD
library(DSOpal)
library(dsBaseClient)
library(dsOmicsClient)
builder <- newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org",
               user = "dsuser", password = "password",
               resource = "test.brge_plink", driver = "OpalDriver")
logindata <- builder$build()

o <- datashield.login(logins = logindata, assign = TRUE,
                          symbol = "client")

ds.ls()
plink.command <- c("--bfile brge --freq")
ans <- ds.PLINK("client", plink.command)
datashield.logout(o)
