package cellgeni.lineagekeep

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.trace.TraceObserverFactoryV2
import nextflow.trace.TraceObserverV2

/**
 * Creates the observer only when `lineageKeep.enabled = true`, so that a config
 * which merely loads the plugin changes nothing.
 */
@CompileStatic
class LineageKeepFactory implements TraceObserverFactoryV2 {

    @Override
    Collection<TraceObserverV2> create(Session session) {
        final opts = session.config.navigate('lineageKeep') as Map ?: Collections.emptyMap()
        final config = new LineageKeepConfig(opts)
        return config.enabled
            ? List.<TraceObserverV2>of(new LineageKeepObserver(config))
            : Collections.<TraceObserverV2>emptyList()
    }
}
