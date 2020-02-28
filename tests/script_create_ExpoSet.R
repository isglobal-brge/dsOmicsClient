rm(list=ls())
library(resourcer)
library(DSOpal)
library(dsBaseClient)
library(DSLite)


library(dsExposome)
library(dsExposomeClient)


# make a DSLite server with resources inside
dslite.server <- newDSLiteServer(resources = list(
  expos = resourcer::newResource(name = "expos", 
                                 url = "https://github.com/isglobal-brge/brgedata/raw/master/data/brge_expo.rda", 
                                 format = "ExposomeSet")
))


# build login details
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "dslite.server", resource = "expos", driver = "DSLiteDriver")
logindata <- builder$build()# login and assign resources
conns <- datashield.login(logins = logindata, assign = TRUE, symbol = "res")# R data file resource


datashield.assign.expr(conns, symbol = "expos", expr = quote(as.resource.object(res)))

ds.class("expos")

ds.familyNames
ds.exposureNames

...


ds.ExWAS


datashield.logout(conns)







builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "test.brgeExpo", driver = "OpalDriver")

logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = TRUE, 
                               symbol = "res")
datashield.assign.expr(conns, symbol = "expo", 
                       expr = quote(as.resource.object(res)))
