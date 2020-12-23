window.onload = function () {
    new Oidc.UserManager().signinRedirectCallback().then(function() {
        window.location = "/"
    })
    
}
