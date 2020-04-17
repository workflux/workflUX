window.onload = function () {
    new Oidc.UserManager({ response_mode: "query"}).signinRedirectCallback().then(function() {
        window.location = "/"
    }).catch(function (e) {
        console.log(e);
        document.getElementById('error_explanation').innerText = "Problem has occured";
        document.getElementById('error_place').innerText = e;
    });   
}
