
			   CEPHEAN-LV2 LIBRARY

cephean-lv2 is a library of LV2 audio plugins developed
in Windows for use in DAWs and the Mod Dwarf effects
pedal. The idea is to make totally rad live sound effects
as fast as I can think them up, so the interface is sparse
and rugged. I'm also trying to make the effects as simple
as possible. Doesn't everyone hate tuning 25 knobs on a
parametric reverb? The Cephean essential verb has 4 knobs,
including wet and dry.

Presently, there are a handful of conventional effects
built in: distortion, reverb, tremelo, etc. and a couple
more notable interesting effects:
	respectrum - fully polyphonic timbre shifter and
		transient shaper with a subband noise gate
	diffuse - transient smearing and blurring,
		generalizing a cabinet/room sim

If you just want a zip file full of plugins, download away!
Otherwise, all the source code and projects are included
here so you can start building your own LV2 plugins in
Windows using my "from scratch" framework built on the LV2
spec itself (no JUCE or DPF overheads). Everything is
included in this single monolithic repo, mostly because I'm
not particularly good at computer science, but ALSO for
convenience!