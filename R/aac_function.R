#' Utils: Area under the curve
#' @description Performs the calculation of the area under the progress curve. Initially created for the plant disease area, whose name is "area under the disease progress curve", it can be adapted to various areas of agrarian science.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @param data Data.frame containing evaluations in columns. Column names must be numeric and not dates or characters
#' @note Just enter the data. Exclude treatment columns. See example.
#' @return Returns a vector with the area values under the curve
#' @references
#'
#' Campbell, C. L., and Madden, L. V. (1990). Introduction to plant disease epidemiology. John Wiley and Sons.
#'
#' @seealso \link{transf}, \link{sketch}
#' @examples
#'
#' #=======================================
#' # Using the simulate1 dataset
#' #=======================================
#' data("simulate1")
#'
#' # Converting to readable format for function
#' dados=cbind(simulate1[simulate1$tempo==1,3],
#'             simulate1[simulate1$tempo==2,3],
#'             simulate1[simulate1$tempo==3,3],
#'             simulate1[simulate1$tempo==4,3],
#'             simulate1[simulate1$tempo==5,3],
#'             simulate1[simulate1$tempo==6,3])
#' colnames(dados)=c(1,2,3,4,5,6)
#' dados
#'
#' # Creating the treatment vector
#' resp=aacp(dados)
#' trat=simulate1$trat[simulate1$tempo==1]
#'
#' # Analyzing by DIC function
#' DIC(trat,resp)
#' @export

aacp=function(data){
  aac <- function(x, y){
    ox <- order(x)
    x <- x[ox]
    y <- y[ox]
    alt <- diff(x)
    bas <- y[-length(y)]+diff(y)/2
    a <- sum(alt*bas)
    return(a)}
  tempo=as.numeric(colnames(data))
  aacp=c(1:nrow(data))
  for(i in 1:nrow(data)){
    aacp[i]=aac(tempo,as.numeric(data[i,]))
  }
  print(aacp)
}
