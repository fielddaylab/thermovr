# VR Thermodynamics Simulator

This is a Unity-based simulation of the thermodynamics properties of water, targeting the Oculus Quest device.

Unity Version: 2019.4.1f1

## Build Steps

* Clone
* Install Correct unity version
*  Android Build Support
* Install Android Studio with
** Android 4.4 (SDK19)
** Android 7.1 (SDK 25)
** SDK Tools (In the “Android SDK Tools” tab)
** SDK - Platform Tools
* Remove any packages that may exist: Legacy XR Input, Oculus (Desktop), OpenVR (Desktop), Oculus XR
* Install "Oculus Integration" asset from App Store. It may ask to upgrade the OVR and Spacailizers and restart Oculus. Do so.
* It will wartn that you need to install the Oculus (Desktop) and OpenVR (Desktop) pacakges. Do so using the package manager.
* In Build Settings, change the platform to Android. MAke sure Scenes/MainScene is selected.

After these steps, the project should build and work, but with the room environment missing.

* Download and Import "3d Living Room" asset from Barking Dog. If it prompts to install package dependancies, do so.
* Do NOT import (or delete if you accidentially imported all) the Barking Dog Scripts wioth the mouse interaction script.
* This will create a bunch of errors. Uninstall the "Package UI" package.


