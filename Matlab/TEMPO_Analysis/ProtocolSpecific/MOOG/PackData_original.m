function [good_data, outputChan] = PackData(good_data,fName)

TEMPO_Defs;

%Load sorted spikes data
load(fName);

%If spsData2 is empty just return
if(isempty(spsData2))
    return;
end

%Find the next channel # to store data in.
chanNumber=4;

%The number of bins and slots per bin
binCount=(spsData2(1).postbuffer+spsData2(1).prebuffer)/spsData2(1).sampleRate*(good_data.htb_header{SPIKE_DB}.speed_units/good_data.htb_header{SPIKE_DB}.speed);
slotsperbin=(spsData2(1).postbuffer+spsData2(1).prebuffer)/binCount;
n=zeros(1,binCount);


for i=1:length(spsData2)
    for j=1:length(spsData2(i).spikeInfo)
        %Bin all the spike times
        edges=[spsData2(i).spikeInfo(j).startCodeTime-spsData2(i).prebuffer+1:slotsperbin:spsData2(i).spikeInfo(j).startCodeTime+spsData2(i).postbuffer+slotsperbin];
        combinedSpikes=spsData2(i).spikeInfo(j).SpikeTimes;
        %combinedSpikes = sort([spsData(index).spikeInfo(j).pSpikeTimes, spsData(index).spikeInfo(j).nSpikeTimes, spsData(index).spikeInfo(j).wSpikeTimes]);
        if(isempty(combinedSpikes))
            n=zeros(1,binCount);
        else
            n=histc(combinedSpikes,edges);
            n=n(1:length(n)-1);
        end
              
        %Put the bins in good_data
        good_data.spike_data(chanNumber+i,:,j)=n;
        outputChan(i)=chanNumber+i;
    end
end

return;