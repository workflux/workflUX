// contains general utilities which are important for multiple modules

function build_loader(message, size="large"){
    if( size == 'large' ){
        loader_class='large_loader'
    } else if( size == 'small' ){
        loader_class='small_loader'
    } else{
        loader_class='tiny_loader'
    }
    loader_html =
        '<div class="w3-container w3-center" style="display:block">' +
                '<div class="' + loader_class + '"></div>' +
                '<div class="w3-center-align">' + message + '</div>' +
        '</div>'
    return loader_html
}


function sleep(secs) {
    secs = (+new Date) + secs * 1000;
    while ((+new Date) < secs);
}