// A little javascript help to propagate values to all empty fields. 

var propagateValues = function (name) {
    var valueToProp = $(`input[name=${name}_propagate]`)[0].value

    if (valueToProp == "") {
        return
    }

    var emptyInputs = $.grep($(`.form-control.${name}`), function (e) { return e.value == "" });
    emptyInputs.forEach(function (i) {
        i.value = valueToProp
    })
}

// genetic, diet, chemical, radiation, temperature, other, control
var toggleGenetic = function () {
    document.querySelectorAll("[id$=_genetic]").forEach(function (e) {
        e.checked = !e.checked;
    });
}

var toggleDiet = function () {
    document.querySelectorAll("[id$=_diet]").forEach(function (e) {
        e.checked = !e.checked;
    });
}

var toggleChemical = function () {
    document.querySelectorAll("[id$=_chemical]").forEach(function (e) {
        e.checked = !e.checked;
    });
}

var toggleRadiation = function () {
    document.querySelectorAll("[id$=_radiation]").forEach(function (e) {
        e.checked = !e.checked;
    });
}

var toggleTemperature = function () {
    document.querySelectorAll("[id$=_temperature]").forEach(function (e) {
        e.checked = !e.checked;
    });
}

var toggleOther = function () {
    document.querySelectorAll("[id$=_other]").forEach(function (e) {
        e.checked = !e.checked;
    });
}

var toggleControl = function () {
    document.querySelectorAll("[id$=_control]").forEach(function (e) {
        e.checked = !e.checked;
    });
}

$(document).ready(function () {
    $("#propButton").click(function () {
        propagateValues("sex")
        propagateValues("devel")
        propagateValues("tissue")
        propagateValues("cell")
        propagateValues("notes")
    });

    $("#toggleGeneticButton").click(function () {
        toggleGenetic()
    });

    $("#toggleDietButton").click(function () {
        toggleDiet()
    });

    $("#toggleChemicalButton").click(function () {
        toggleChemical()
    });

    $("#toggleRadiationButton").click(function () {
        toggleRadiation()
    });

    $("#toggleTemperatureButton").click(function () {
        toggleTemperature()
    });

    $("#toggleOtherButton").click(function () {
        toggleOther()
    });

    $("#toggleControlButton").click(function () {
        toggleControl()
    });

});

