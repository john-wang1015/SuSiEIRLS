#' @export
main_tptn=function(true_main_index,main_index){

# if there is no true main effect, tp = tn = 1 only if main_index is also null
if(is.null(true_main_index)){
tp=tn=ifelse(is.null(main_index)==1,1,0)
}

tp=ifelse(length(setdiff(true_main_index,main_index$Index))==0,1,0)
true_main_cs=unique(main_index$CS[match(true_main_index,main_index$Index)])
tn = ifelse(length(setdiff(main_index$CS, true_main_cs)) == 0, 1, 0)

return(list(tp=tp,tn=tn))
}
