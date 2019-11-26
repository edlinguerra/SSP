\name{epibionts}
\alias{epibionts}
\docType{data}
\title{
Epibionts on Caribbean mangrove roots
}
\description{
Data corresponds to epibenthic organisms in mangrove roots from Laguna de La Restinga National Park, Venezuela (Guerra-Castro et al. 2016). Data consists of the coverage (by point-intercept) of 110 taxa identified in 240 mangrove roots, sampled under a hierarchically nested spatial design that included four random sites within each of three sectors of the lagoon system corresponding to a strong environmental gradient: external (E), intermediate (M), and internal (I). The abundance of epibenthic organisms of 8 roots were described within each site, producing a total of 32 roots in each sector. This spatial protocol was repeated five times over a period of 14 months. For demonstrative purpose, data from the 4th sampling period was randomly chosen as data for this package.
}
\usage{data("epibionts")}
\format{
  A data frame with 96 observations on the following 152 variables.
  \describe{
    \item{\code{sector}}{a factor with levels \code{E} \code{I} \code{M}}
    \item{\code{site}}{a numeric vector}
    \item{\code{Aaptos.sp}}{a numeric vector}
    \item{\code{Acanthophora.spicifera}}{a numeric vector}
    \item{\code{Acetabularia.crenulata}}{a numeric vector}
    \item{\code{Aglaothamnion.sp}}{a numeric vector}
    \item{\code{Amathia.sp}}{a numeric vector}
    \item{\code{Amorphinopsis.atlantica}}{a numeric vector}
    \item{\code{Amphimedon.erina}}{a numeric vector}
    \item{\code{Anemonia.sargassensis}}{a numeric vector}
    \item{\code{Aplidium.accarense}}{a numeric vector}
    \item{\code{Aplysilla.glacialis}}{a numeric vector}
    \item{\code{Ascidia.curvata}}{a numeric vector}
    \item{\code{Ascidia.sp}}{a numeric vector}
    \item{\code{Ascidia.sydneiensis}}{a numeric vector}
    \item{\code{Balanus.sp}}{a numeric vector}
    \item{\code{Bartholomea.annulata}}{a numeric vector}
    \item{\code{Biemna.caribea}}{a numeric vector}
    \item{\code{Bostrychia.tenella}}{a numeric vector}
    \item{\code{Botrylloides.nigrum}}{a numeric vector}
    \item{\code{Botrylloides.sp.1}}{a numeric vector}
    \item{\code{Botrylloides.sp.2}}{a numeric vector}
    \item{\code{Brachidontes.exustus}}{a numeric vector}
    \item{\code{Branchiomma.conspersum}}{a numeric vector}
    \item{\code{Branchiomma.nigromaculatum}}{a numeric vector}
    \item{\code{Bryopsis.sp}}{a numeric vector}
    \item{\code{Bugula.neritina}}{a numeric vector}
    \item{\code{Bugula.sp}}{a numeric vector}
    \item{\code{Calliactis.tricolor}}{a numeric vector}
    \item{\code{Callyspongia..Callyspongia..pallida}}{a numeric vector}
    \item{\code{Carijoa.riisei}}{a numeric vector}
    \item{\code{Caulerpa.racemosa}}{a numeric vector}
    \item{\code{Caulerpa.racemosa.var.peltata}}{a numeric vector}
    \item{\code{Caulerpa.sertularioides}}{a numeric vector}
    \item{\code{Caulerpa.verticillata}}{a numeric vector}
    \item{\code{Caulibugula.sp}}{a numeric vector}
    \item{\code{Celleporaria.sp}}{a numeric vector}
    \item{\code{Ceramium.diaphanum}}{a numeric vector}
    \item{\code{Chaetomorpha.sp.1}}{a numeric vector}
    \item{\code{Chaetomorpha.sp.2}}{a numeric vector}
    \item{\code{Chalinula.molitba}}{a numeric vector}
    \item{\code{Chelonaplysilla.erecta}}{a numeric vector}
    \item{\code{Chondrilla.nucula}}{a numeric vector}
    \item{\code{Chthamalus.sp}}{a numeric vector}
    \item{\code{Clathria..Clathria..microchela}}{a numeric vector}
    \item{\code{Clathria.sp}}{a numeric vector}
    \item{\code{Clavelina.oblonga}}{a numeric vector}
    \item{\code{Clavelina.picta}}{a numeric vector}
    \item{\code{Complejo.Cliona.celata}}{a numeric vector}
    \item{\code{Crassostrea.rhizophorae}}{a numeric vector}
    \item{\code{Dictyota.sp}}{a numeric vector}
    \item{\code{Didemnum.cineraceum}}{a numeric vector}
    \item{\code{Didemnum.perlucidum}}{a numeric vector}
    \item{\code{Didemnum.sp}}{a numeric vector}
    \item{\code{Diplosoma.listerianum}}{a numeric vector}
    \item{\code{Distaplia.bermudensis}}{a numeric vector}
    \item{\code{Distaplia.stylifera}}{a numeric vector}
    \item{\code{Dynamena.sp}}{a numeric vector}
    \item{\code{Dysidea.etheria}}{a numeric vector}
    \item{\code{Dysidea.sp}}{a numeric vector}
    \item{\code{Ecteinascidia.sp}}{a numeric vector}
    \item{\code{Ecteinascidia.styeloides}}{a numeric vector}
    \item{\code{Ecteinascidia.turbinata}}{a numeric vector}
    \item{\code{Eudistoma.olivaceum}}{a numeric vector}
    \item{\code{Eusynstyela.tincta}}{a numeric vector}
    \item{\code{Exaiptasia.pallida}}{a numeric vector}
    \item{\code{Ficopomatus.sp}}{a numeric vector}
    \item{\code{Geodia.papyracea}}{a numeric vector}
    \item{\code{Halichondria..Halichondria..magniconulosa}}{a numeric vector}
    \item{\code{Halichondria..Halichondria..melanadocia}}{a numeric vector}
    \item{\code{Haliclona..Halichoclona..magnifica}}{a numeric vector}
    \item{\code{Haliclona..Reniera..implexiformis}}{a numeric vector}
    \item{\code{Haliclona..Reniera..manglaris}}{a numeric vector}
    \item{\code{Haliclona..Reniera..ruetzleri}}{a numeric vector}
    \item{\code{Haliclona..Reniera..tubifera}}{a numeric vector}
    \item{\code{Haliclona..Rhizoniera..curacaoensis}}{a numeric vector}
    \item{\code{Haliclona..Soestella..caerulea}}{a numeric vector}
    \item{\code{Haliclona..Soestella..smithae}}{a numeric vector}
    \item{\code{Haliclona..Soestella..twincayensis}}{a numeric vector}
    \item{\code{Halimeda.sp}}{a numeric vector}
    \item{\code{Halisarca.sp}}{a numeric vector}
    \item{\code{Halopteris.sp}}{a numeric vector}
    \item{\code{Herdmania.pallida}}{a numeric vector}
    \item{\code{Hippopodina.feegeensis}}{a numeric vector}
    \item{\code{Hydroides.sp}}{a numeric vector}
    \item{\code{Hyrtios.proteus}}{a numeric vector}
    \item{\code{Iotrochota.birotulata}}{a numeric vector}
    \item{\code{Ircinia.felix}}{a numeric vector}
    \item{\code{Ircinia.sp}}{a numeric vector}
    \item{\code{Isognomon.alatus}}{a numeric vector}
    \item{\code{Kirchenpaueria.sp}}{a numeric vector}
    \item{\code{Lissoclinum.sp}}{a numeric vector}
    \item{\code{Lissodendoryx..Lissodendoryx..isodictyalis}}{a numeric vector}
    \item{\code{Lithophyllum.pustulatum}}{a numeric vector}
    \item{\code{Microcosmus.exasperatus}}{a numeric vector}
    \item{\code{Molgula.occidentalis}}{a numeric vector}
    \item{\code{Murrayella.periclados}}{a numeric vector}
    \item{\code{Mycale..Aegogropila..carmigropila}}{a numeric vector}
    \item{\code{Mycale..Aegogropila..citrina}}{a numeric vector}
    \item{\code{Mycale..Carmia..magnirhaphidifera}}{a numeric vector}
    \item{\code{Mycale..Carmia..microsigmatosa}}{a numeric vector}
    \item{\code{Mycale..Mycale..laevis}}{a numeric vector}
    \item{\code{Mycale..Zygomycale..angulosa}}{a numeric vector}
    \item{\code{Mycale.sp}}{a numeric vector}
    \item{\code{Nemalecium.sp}}{a numeric vector}
    \item{\code{Notaulax.nudicollis}}{a numeric vector}
    \item{\code{Obelia.sp}}{a numeric vector}
    \item{\code{Oceanapia.nodosa}}{a numeric vector}
    \item{\code{Padina.sp}}{a numeric vector}
    \item{\code{Perna.viridis}}{a numeric vector}
    \item{\code{Perophora.viridis}}{a numeric vector}
    \item{\code{Phaeophyceae}}{a numeric vector}
    \item{\code{Phallusia.nigra}}{a numeric vector}
    \item{\code{Phyllangia.americana}}{a numeric vector}
    \item{\code{Pinctada.imbricata}}{a numeric vector}
    \item{\code{Plakortis.angulospiculatus}}{a numeric vector}
    \item{\code{Polyclinum.constellatum}}{a numeric vector}
    \item{\code{Polysiphonia.sp.1}}{a numeric vector}
    \item{\code{Polysiphonia.sp.3}}{a numeric vector}
    \item{\code{Polysiphonia.subtilissima}}{a numeric vector}
    \item{\code{Pteria.colymbus}}{a numeric vector}
    \item{\code{Pyura.sp..1}}{a numeric vector}
    \item{\code{Pyura.sp..2}}{a numeric vector}
    \item{\code{Pyura.vittata}}{a numeric vector}
    \item{\code{Rhizoclonium.sp}}{a numeric vector}
    \item{\code{Rhodosoma.turcicum}}{a numeric vector}
    \item{\code{Sabella.sp}}{a numeric vector}
    \item{\code{Sabellastarte.magnifica}}{a numeric vector}
    \item{\code{Schizoporella.pungens}}{a numeric vector}
    \item{\code{Scopalina.ruetzleri}}{a numeric vector}
    \item{\code{Scopalina.sp}}{a numeric vector}
    \item{\code{Scrupocellaria.sp}}{a numeric vector}
    \item{\code{Sphacelaria.rigidula}}{a numeric vector}
    \item{\code{Spongia..Spongia..pertusa}}{a numeric vector}
    \item{\code{Spongia..Spongia..tubulifera}}{a numeric vector}
    \item{\code{Sporolithon.episporum}}{a numeric vector}
    \item{\code{Spyridia.hypnoides}}{a numeric vector}
    \item{\code{Styela.canopus}}{a numeric vector}
    \item{\code{Styela.sp.1}}{a numeric vector}
    \item{\code{Styela.sp.2}}{a numeric vector}
    \item{\code{Suberites.aurantiacus}}{a numeric vector}
    \item{\code{Symplegma.brakenhielmi}}{a numeric vector}
    \item{\code{Symplegma.rubra}}{a numeric vector}
    \item{\code{Synnotum.circinatum}}{a numeric vector}
    \item{\code{Tedania..Tedania..ignis}}{a numeric vector}
    \item{\code{Terpios.manglaris}}{a numeric vector}
    \item{\code{Tethya.actinia}}{a numeric vector}
    \item{\code{Tethya.sp}}{a numeric vector}
    \item{\code{Trididemnum.orbiculatum}}{a numeric vector}
    \item{\code{Ulva.sp}}{a numeric vector}
    \item{\code{Viatrix.globulifera}}{a numeric vector}
    \item{\code{Zoobotryon.verticillatum}}{a numeric vector}
  }
}
\details{
Data corresponds to epibenthic organisms in mangrove roots from Laguna de La Restinga National Park, Venezuela (Guerra-Castro et al. 2016). Data consists of the coverage (by point-intercept) of 110 taxa identified in 240 mangrove roots, sampled under a hierarchically nested spatial design that included four random sites within each of three sectors of the lagoon system corresponding to a strong environmental gradient: external (E), intermediate (M), and internal (I). The abundance of epibenthic organisms of 8 roots were described within each site, producing a total of 32 roots in each sector. This spatial protocol was repeated five times over a period of 14 months. For demonstrative purpose, data from the 4th sampling period was randomly chosen as data for this package.
}
\source{
https://doi.org/10.3354/meps11693
}
\references{
Guerra-Castro, E. J., J. E. Conde, and J. J. Cruz-Motta. (2016). Scales of spatial variation in tropical benthic assemblages and their ecological relevance: epibionts on Caribbean mangrove roots as a model system. Marine Ecology Progress Series 548:97-110.
}
\examples{
data(epibionts)
str(epibionts)
}
\keyword{datasets}
\keyword{SSP}
