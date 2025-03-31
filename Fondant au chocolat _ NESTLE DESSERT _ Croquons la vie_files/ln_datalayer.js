/**
 * @file
 * Initializes some client side datalayer tags.
 */
function getCookie(cname) {
  'use strict';
  var name = cname + '=';
  var ca = document.cookie.split(';');
  for (var i = 0; i < ca.length; i++) {
    var c = ca[i];
    while (c.charAt(0) === ' ') {
      c = c.substring(1);
    }
    if (c.indexOf(name) === 0) {
      return c.substring(name.length, c.length);
    }
  }
  return '';
}

if (typeof dataLayer_tags !== 'undefined' && dataLayer_tags.userInformation) {
  if (typeof dataLayer_tags.userInformation.gaClientID !== 'undefined') {
    dataLayer_tags.userInformation.gaClientID = getCookie('_ga');
  }
  if (typeof dataLayer_tags.userInformation.deviceType !== 'undefined') {
    if (typeof isMobile !== undefined) {
      var deviceType = 'Desktop';
      if (isMobile.tablet) {
        deviceType = 'Tablet';
      }
      if (isMobile.phone) {
        deviceType = 'Mobile';
      }
      dataLayer_tags.userInformation.deviceType = deviceType;
    }
  }
}
