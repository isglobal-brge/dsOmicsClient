library(DSOpal)
library(dsBaseClient)
library(dsOmics)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", user = "dsuser", password = "password", resource = "test.gtex_kidney", driver = "OpalDriver")
logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "res")

ds.ls(conns)


datashield.assign.expr(conns, symbol = "ES", expr = quote(as.resource.object(res)))



library(opalr)
o <- opal.login(username="administrator", password="password", url="https://opal-test.obiba.org")
opal.execute(o, "BiocManager::install('SummarizedExperiment', ask=FALSE)")


library(opalr)
o <- opal.login(username="administrator", password="password", url="https://opal-test.obiba.org")
opal.execute(o, "BiocManager::install('MEAL', ask=FALSE)")
