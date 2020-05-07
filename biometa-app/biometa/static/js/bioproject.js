// A little javascript help to propagate values to all empty fields. 

var propagateValues = function(name){
    var valueToProp = $(`input[name=${name}_propagate]`)[0].value

    if (valueToProp == "") {
        return
    }

    var emptyInputs = $.grep($(`.form-control.${name}`), function(e){return e.value == ""});
    emptyInputs.forEach(function(i) {
        i.value = valueToProp
    })
}

$(document).ready(function(){
    $("#propButton").click(function(){
        propagateValues("sex")
        propagateValues("dev")
        propagateValues("tissue")
        propagateValues("cell")
    });
});

